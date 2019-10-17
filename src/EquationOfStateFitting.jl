"""
# module EquationOfStateFitting



# Examples

```jldoctest
julia>
```
"""
module EquationOfStateFitting

using LinearAlgebra: det

using Dates: Second
using Compat: isnothing
using EquationsOfState
using EquationsOfState.Collections: EquationOfState
using EquationsOfState.NonlinearFitting: lsqfit
using EquationsOfState.Find: findvolume
using ExtensibleScheduler
using Kaleido: @batchlens
using MLStyle: @match
using QuantumESPRESSOBase: to_qe, option
using QuantumESPRESSOBase.Inputs.PWscf: PWscfInput, autofill_cell_parameters
using QuantumESPRESSOParsers.OutputParsers.PWscf: parse_total_energy,
                                                  parse_cell_parameters,
                                                  parse_head,
                                                  isjobdone
using Setfield: set
using Unitful
using Unitful: AbstractQuantity
using UnitfulAtomic

import ..Step
using ..SelfConsistentField: write_metadata

export update_alat_press, prepare, finish, submit, query, write_job, isdone, checkjob

function update_alat_press(
    template::PWscfInput,
    eos::EquationOfState{<:AbstractQuantity},
    pressure::AbstractQuantity,
)
    if isnothing(template.cell_parameters)
        template = autofill_cell_parameters(template)
    end
    # In case `eos.v0` has a `Int` as `T`. See https://github.com/PainterQubits/Unitful.jl/issues/274.
    v0 = float(eos.v0)
    volume = findvolume(PressureForm(), eos, pressure, (eps(v0), 1.3v0))
    # If the `CellParametersCard` contains a matrix of plain numbers (no unit).
    determinant = @match option(template.cell_parameters) begin
        # `alat` uses relative values WRT `celldm`, which uses "bohr" as unit.
        # So `"alat"` is equivalent to `"bohr"`.
        "alat" || "bohr" => det(template.cell_parameters.data) * u"bohr^3"
        "angstrom" => det(template.cell_parameters.data) * u"angstrom^3" |> u"bohr^3"
    end
    # `cbrt` works with units.
    alat = cbrt(volume / determinant) |> NoUnits  # This is dimensionless.
    lenses = @batchlens(begin
        _.system.celldm  # Get the `template`'s `system.celldm` value
        _.cell.press     # Get the `template`'s `cell.press` value
        _.cell_parameters.option
    end)
    # Set the `template`'s `system.celldm` and `cell.press` values
    return set(template, lenses, ([alat], ustrip(u"kbar", pressure), "alat"))
end # function update_alat_press

# This is a helper function and should not be exported.
function _boilerplate(step::Step{N}, template::PWscfInput) where {N}
    lenses = @batchlens(begin
        _.control.calculation  # Get the `template`'s `control.calculation` value
        _.control.verbosity    # Get the `template`'s `control.verbosity` value
        _.control.tstress      # Get the `template`'s `control.tstress` value
        _.control.tprnfor      # Get the `template`'s `control.tprnfor` value
    end)
    # Set the `template`'s values with...
    return set(template, lenses, (N == 1 ? "scf" : "vc-relax", "high", true, true))
end # function _boilerplate

function prepare(
    step::Step,
    inputs::AbstractVector{<:AbstractString},
    template::PWscfInput,
    trial_eos::EquationOfState,
    pressures::AbstractVector,
    metadatafiles::AbstractVector{<:AbstractString} = map(x -> splitext(x)[1] * ".json", inputs),
    verbose::Bool = false
)
    # Check parameters
    @assert(
        length(inputs) == length(pressures) == length(metadatafiles),
        "The inputs, pressures and the metadata files must be the same size!"
    )
    template = _boilerplate(step, template)
    # Write input and metadata
    for (input, pressure, metadatafile) in zip(inputs, pressures, metadatafiles)
        # Get a new `object` from the `template`, with its `alat` and `pressure` changed
        object = update_alat_press(template, trial_eos, pressure)
        write(input, to_qe(object, verbose = verbose))  # Write the `object` to a Quantum ESPRESSO input file
        write_metadata(metadatafile, input, template)
    end
    return
end # function prepare

function finish(
    ::Step{N},
    outputs::AbstractVector{<:AbstractString},
    trial_eos::EquationOfState{<:AbstractQuantity},
) where {N}
    energies, volumes = zeros(length(outputs)), zeros(length(outputs))
    for (i, output) in enumerate(outputs)
        open(output, "r") do io
            s = read(io, String)
            isjobdone(s) || @warn("Job is not finished!")
            energies[i] = parse_total_energy(s)[end]
            volumes[i] = @match N begin
                1 => parse_head(s)["unit-cell volume"]
                2 => det(parse_cell_parameters(s)[end])
                _ => error("The step $N must be `1` or `2`!")
            end
        end
    end
    return lsqfit(EnergyForm(), trial_eos, volumes .* u"bohr^3", energies .* u"Ry")
end # function finish

mutable struct CalculationStatus{T}
    io::AbstractDict{T,T}
    pending::AbstractDict{T,T}
    running::AbstractDict{T,T}
    done::AbstractDict{T,T}
end

function communicate(cmd::Cmd)
    out = Pipe()
    err = Pipe()

    process = run(pipeline(cmd, stdout = out, stderr = err), wait = false)
    close(out.in)
    close(err.in)

    stdout = @async String(read(out))
    stderr = @async String(read(err))
    wait(process)
    return (
        stdout = fetch(stdout),
        stderr = fetch(stderr),
        code = process.exitcode
    )
end

function write_job(filename::AbstractString, io::AbstractDict{T,T}, account::AbstractString) where {T<:AbstractString}
    n = div(24, length(io))
    array = join(["[\"$k\"]=\"$v\"" for (k, v) in io], " ")
    s = """
    #!/usr/bin/env bash

    #SBATCH -A $account
    #SBATCH -N 1
    #SBATCH --tasks-per-node=24
    #SBATCH -J job
    #SBATCH --time=01:00:00

    module load intel-parallel-studio/2017

    declare -A array=($array)

    for i in "\${!array[@]}"; do
        mpirun -np $n pw.x -in "\$i" > "\${array[\$i]}" &
        sleep 3
    done
    wait
    """
    write(filename, s)
    return
end # function write_job

function submit(job::T) where {T<:AbstractString}
    out, err, code = communicate(`sbatch --parsable $job > jobid`)
    if code == 0
        return chomp(out)  # It is the job id.
    else
        error(err)
    end
end # function submit

function query(jobid::AbstractString)
    out, err, code = communicate(`squeue --job=$jobid`)
    if code != 0
        error(err)
    end
    for line in split(out, '\n')
        m = match(Regex(jobid), line)
        !isnothing(m) && return strip(line)
    end
    return
end # function query

function isdone(jobid::AbstractString)
    return isnothing(query(jobid)) ? true : false
end # function isdone

function checkjob(jobid::AbstractString, dt::Int = 10)
    # From https://scls19fr.github.io/ExtensibleScheduler.jl/dev/getting_started/
    function f(x, sch)
        isdone(x) && ExtensibleScheduler.shutdown(sch)
    end
    sched = BlockingScheduler()
    action = Action(f, jobid, sched)
    trigger = Trigger(Second(dt), n = -1)  # Infinte loop until job is done
    add(sched, action, trigger)
    run(sched)
end # function overall

# function workflow(
#     io::AbstractDict{T,T},
#     template::PWscfInput,
#     trial_eos::EquationOfState,
#     pressures::AbstractVector,
#     metadatafiles::AbstractVector{<:AbstractString} = map(x -> splitext(x)[1] * ".json", keys(io)),
#     account::AbstractString,
#     verbose::Bool = false
# ) where {T<:AbstractString}
#     prepare(Step(1), keys(io), template, trial_eos, pressures, metadatafiles)
#     write_job("job.sh", io, account)
#     jobid = submit("job.sh")
#     isdone(jobid, account)
#     trial_eos = finish(values(io), trial_eos)
#     # prepare(Step(2), keys(io), template, trial_eos, pressures, metadatafiles)
#     # write_job("job.sh", io, account)
#     # jobid = submit("job.sh")
#     # return finish(values(io), trial_eos)
# end # function workflow

end
