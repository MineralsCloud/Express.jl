"""
# module StructureOptimization



# Examples

```jldoctest
julia>
```
"""
module StructureOptimization

using LinearAlgebra

using EquationsOfState
using EquationsOfState.Collections
using EquationsOfState.NonlinearFitting
using EquationsOfState.FindVolume
using IntervalArithmetic
using IntervalRootFinding
import JSON
using Kaleido: @batchlens
using QuantumESPRESSOBase
using QuantumESPRESSOBase.Inputs.PWscf
using QuantumESPRESSOParsers.InputParsers
using Setfield
using SlurmWorkloadFileGenerator.Commands
using SlurmWorkloadFileGenerator.SystemModules
using SlurmWorkloadFileGenerator.Scriptify
using SlurmWorkloadFileGenerator.Shells

using Express

export update_alat_press, write_input, write_script, write_metadata, prepare, finish

function update_alat_press(template::PWscfInput, eos::EquationOfState, pressure::Real)
    volume = find_volume(PressureRelation, eos, pressure, 0..1000, Newton).interval.lo
    alat = cbrt(volume / det(template.cell_parameters.data))
    lenses = @batchlens begin
        _.system.celldm ∘ _[$1]
        _.cell.press
    end
    set(template, lenses, (alat, pressure))
end # function update_alat_press

function write_input(
    input::AbstractString,
    template::PWscfInput,
    eos::EquationOfState,
    pressure::Real,
    verbose::Bool = false
)
    object = update_alat_press(template, eos, pressure)
    open(input, "r+") do io
        write(io, to_qe(object, verbose = verbose))
    end
end # function write_input
function write_input(
    inputs::AbstractVector{<:AbstractString},
    template::PWscfInput,
    eos::EquationOfState,
    pressures::AbstractVector{<:Real},
    verbose::Bool = false
)
    length(inputs) == length(pressures) || throw(DimensionMismatch("The number of inputs should equal the number of pressures!"))
    [write_input(input, template, eos, pressure, verbose) for (input, pressure) in zip(inputs, pressures)]
end # function write_input

function write_script(shell::Shell, sbatch::Sbatch, modules, pressures::AbstractVecOrMat)
    content = "#!$(string(shell.path))\n"
    content *= scriptify(sbatch) * '\n'
    content *= scriptify(modules) * '\n'
    content *= """
    for i in $(join(pressures, ' ')); do
        (
            cd vc_\$i || return
            mpirun -np $(div(24, length(pressures))) pw.x -in vc_\$i.in > vc_\$i.out &
            sleep 3
        )
    done
    wait
    """
    open(sbatch.script, "r+") do io
        write(io, content)
    end
end # function write_script

function write_metadata(output::AbstractString, object::PWscfInput, input::AbstractString)
    metadata = Dict(
        "outdir" => object.control.outdir,
        "prefix" => object.control.prefix,
        "pseudo_dir" => object.control.pseudo_dir,
        "pseudopotentials" => [getfield(x, :pseudopotential) for x in object.atomic_species.data],
        "input" => input
    )
    metadata["wfcdir"] = if object.control.wf_collect
        metadata["outdir"]
    else
        object.control.wfcdir
    end
    if object.control.lkpoint_dir
        metadata["lkpoint_dir"] = metadata["outdir"] * "/" * metadata["prefix"] * ".save"
    end
    metadata["wfc_namepattern"] = metadata["prefix"] * ".wfc"
    lowercase(splitext(output)[2]) != ".json" && error("The file to be dumped must be a JSON file!")
    open(output, "r+") do io
        JSON.print(io, metadata)
    end
end # function write_metadata

prepare(::Step, args...) = error("The only allowed steps are `1` and `2`!")
function prepare(
    step::Step{1},
    inputs::AbstractVector,
    template::PWscfInput,
    trial_eos::EquationOfState,
    pressures::AbstractVector,
    metadatafiles::AbstractVector
)
    isnothing(template.cell_parameters) && (template = autofill_cell_parameters(template))
    if template.control.calculation != "scf"
        @warn "The calculation type is $(template.control.calculation), not \"scf\"! We will set it for you."
        @set! template.control.calculation = "scf"
    end
    write_input(inputs, template, trial_eos, pressures)
    write_script(
        Shell("/bin/bash"),
        Sbatch([NTASKS_PER_NODE(24), TIME("24:00:00"), NODES(2)], "job.sh"),
        [SystemModule("intel-parallel-studio/2017")],
        pressures
    )
    write_metadata.(metadatafiles, template, inputs)
end # function prepare
function prepare(
    step::Step{2},
    new_inputs::AbstractVector,
    previous_outputs::AbstractVector,
    template::PWscfInput,
    trial_eos::EquationOfState,
    pressures::AbstractVector,
    volumes::AbstractVector,
    metadatafiles::AbstractVector
)
    length(new_inputs) == length(previous_outputs) == length(pressures) ==
    length(volumes) && throw(DimensionMismatch("The previous inputs, new inputs, the pressures to be applied, and the volumes of that must have the same length!"))
    if template.control.calculation != "vc-relax"
        @warn "The calculation type is $(template.control.calculation), not \"vc-relax\"! We will set it for you."
        @set! template.control.calculation = "vc-relax"
    end
    isnothing(template.cell_parameters) && (template = autofill_cell_parameters(template))
    energies = parse_total_energy.(previous_outputs)
    eos = lsqfit(EnergyTarget, trial_eos, volumes, energies)
    write_input(new_inputs, template, eos, pressures)
    write_script(
        Shell("/bin/bash"),
        Sbatch([NTASKS_PER_NODE(24), TIME("24:00:00"), NODES(2)], "job.sh"),
        [SystemModule("intel-parallel-studio/2017")],
        pressures
    )
    write_metadata.(metadatafiles, template, new_inputs)
end # function prepare

function finish(outputs::AbstractVector, trial_eos, volumes::AbstractVector)
    energies = parse_total_energy.(outputs)
    eos = lsqfit(EnergyTarget, trial_eos, volumes, energies)
    return eos
end # function finish

end
