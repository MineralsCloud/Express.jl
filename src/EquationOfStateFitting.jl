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
using QuantumESPRESSO: to_qe, cell_volume
using QuantumESPRESSO.Cards.PWscf: AtomicPositionsCard, CellParametersCard
using QuantumESPRESSO.Inputs: autofill_cell_parameters
using QuantumESPRESSO.Inputs.PWscf: PWInput
using QuantumESPRESSO.Outputs.PWscf: Preamble,
                                     parse_electrons_energies,
                                     parsefinal,
                                     isjobdone
using Setfield: set
using Unitful
using UnitfulAtomic

import ..Step
using ..SelfConsistentField: write_metadata

export update_alat_press, prepare, finish, submit

function update_alat_press(
    template::PWInput,
    eos::EquationOfState{<:Unitful.AbstractQuantity},
    pressure::Unitful.AbstractQuantity,
)
    if isnothing(template.cell_parameters)
        template = autofill_cell_parameters(template)
    end
    # In case `eos.v0` has a `Int` as `T`. See https://github.com/PainterQubits/Unitful.jl/issues/274.
    v0 = float(eos.v0)
    volume = findvolume(PressureForm(), eos, pressure, (eps(v0), 1.3v0))
    determinant = cell_volume(template) * u"bohr^3"
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
function _boilerplate(step::Step{N}, template::PWInput) where {N}
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
    template::PWInput,
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
    trial_eos::EquationOfState{<:Unitful.AbstractQuantity},
) where {N}
    energies, volumes = zeros(length(outputs)), zeros(length(outputs))
    for (i, output) in enumerate(outputs)
        open(output, "r") do io
            s = read(io, String)
            isjobdone(s) || @warn("Job is not finished!")
            energies[i] = parse_electrons_energies(s, :combined)[end]
            volumes[i]= if N == 1
                parse(Preamble, s)["unit-cell volume"]
            elseif N == 2
                det(parsefinal(CellParametersCard, s)[end])
            else
                error("The step $N must be `1` or `2`!")
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

function submit(job::T) where {T<:AbstractString}
    out, err, code = communicate(`sbatch --parsable $job > jobid`)
    if code == 0
        return chomp(out)  # It is the job id.
    else
        error(err)
    end
end # function submit

# function workflow(
#     io::AbstractDict{T,T},
#     template::PWInput,
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
