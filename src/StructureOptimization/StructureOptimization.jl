"""
# module StructureOptimization



# Examples

```jldoctest
julia>
```
"""
module StructureOptimization

using LinearAlgebra

import JSON
using Kaleido: @batchlens

using IntervalArithmetic
using IntervalRootFinding
using EquationsOfState
using EquationsOfState.Collections
using EquationsOfState.NonlinearFitting
using EquationsOfState.NumericallyFindVolume
using QuantumESPRESSOBase
using QuantumESPRESSOBase.Inputs.PWscf
using QuantumESPRESSOParsers.InputParsers
using Setfield
using SlurmWorkloadFileGenerator.Commands
using SlurmWorkloadFileGenerator.SystemModules
using SlurmWorkloadFileGenerator.Scriptify
using SlurmWorkloadFileGenerator.Shells

export update_alat_press, generate_input!, generate_script, dump_metadata!, prepare_for_step, postprocess

function update_alat_press(template::PWscfInput, eos::EquationOfState, pressure::Real)
    volume = find_volume(PressureTarget, eos, pressure, 0..1000, Newton).interval.lo
    alat = cbrt(volume / det(template.cell_parameters.data))
    lenses = @batchlens begin
        _.system.celldm ∘ _[$1]
        _.cell.press
    end
    set(template, lenses, (alat, pressure))
end # function update_alat_press

function generate_input!(
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
end # function generate_input!
function generate_input!(
    inputs::AbstractVector{<:AbstractString},
    template::PWscfInput,
    eos::EquationOfState,
    pressures::AbstractVector{<:Real},
    verbose::Bool = false
)
    length(inputs) == length(pressures) || throw(DimensionMismatch("The number of inputs should equal the number of pressures!"))
    [generate_input!(input, template, eos, pressure, verbose) for (input, pressure) in zip(inputs, pressures)]
end # function generate_input!

function generate_script(shell::Shell, sbatch::Sbatch, modules, pressures::AbstractVecOrMat)
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
end # function generate_script

function dump_metadata!(output::AbstractString, object::PWscfInput, input::AbstractString)
    metadata = Dict(
        "outdir" => object.control.outdir,
        "prefix" => object.control.prefix,
        "lkpoint_dir" => object.control.lkpoint_dir,
        "pseudo_dir" => object.control.pseudo_dir,
        "pseudopotentials" => map(getfield(:pseudopotential), object.atomic_species.data),
        "input" => input
    )
    if object.control.wf_collect
        metadata["wfcdir"] = object.control.wfcdir
    end
    splitext(output) != "json" && error("The file to be dumped must be a JSON file!")
    open(output, "r+") do io
        JSON.print(io, metadata)
    end
end # function dump_metadata!

prepare_for_step(step::Int, args...) = prepare_for_step(Val(step), args...)
function prepare_for_step(
    step::Val{1},
    inputs::AbstractVector,
    template::PWscfInput,
    trial_eos::EquationOfState,
    pressures::AbstractVector,
    metadatafiles::AbstractVector
)
    isnothing(template.cell_parameters) && (template = autogenerate_cell_parameters(template))
    generate_input!(inputs, template, trial_eos, pressures)
    generate_script(
        Shell("/bin/bash"),
        Sbatch([NTASKS_PER_NODE(24), TIME("24:00:00"), NODES(2)], "job.sh"),
        [SystemModule("intel-parallel-studio/2017")],
        pressures
    )
    dump_metadata!.(metadatafiles, template, inputs)
end # function prepare_for_step
function prepare_for_step(
    step::Val{2},
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
    isnothing(template.cell_parameters) && (template = autogenerate_cell_parameters(template))
    energies = parse_total_energy.(previous_outputs)
    eos = lsqfit(EnergyTarget, trial_eos, volumes, energies)
    generate_input!(new_inputs, template, eos, pressures)
    generate_script(
        Shell("/bin/bash"),
        Sbatch([NTASKS_PER_NODE(24), TIME("24:00:00"), NODES(2)], "job.sh"),
        [SystemModule("intel-parallel-studio/2017")],
        pressures
    )
    dump_metadata!.(metadatafiles, template, new_inputs)
end # function prepare_for_step

function postprocess(outputs::AbstractVector, trial_eos, volumes::AbstractVector)
    energies = parse_total_energy.(outputs)
    eos = lsqfit(EnergyTarget, trial_eos, volumes, energies)
    return eos
end # function postprocess

end
