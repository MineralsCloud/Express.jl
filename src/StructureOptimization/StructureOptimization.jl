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

export update_alat_press, write_input, write_metadata, prepare, finish

function update_alat_press(template::PWscfInput, eos::EquationOfState, pressure::Real)
    volume = find_volume(PressureForm(), eos, pressure, 0..1000, Newton).interval.lo
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

function write_metadata(output::AbstractString, object::PWscfInput, input::AbstractString)
    metadata = Dict(
        "outdir" => object.control.outdir,
        "prefix" => object.control.prefix,
        "pseudo_dir" => object.control.pseudo_dir,
        "pseudopotentials" => [getfield(x, :pseudopotential) for x in object.atomic_species.data],
        "input" => input
    )
    metadata["wfcdir"] = object.control.wf_collect ? metadata["outdir"] : object.control.wfcdir
    metadata["wfc_namepattern"] = metadata["prefix"] * ".wfc"
    if object.control.lkpoint_dir
        metadata["lkpoint_dir"] = metadata["outdir"] * "/" * metadata["prefix"] * ".save"
    end
    lowercase(splitext(output)[2]) != ".json" && error("The file to be dumped must be a JSON file!")
    open(output, "r+") do io
        JSON.print(io, metadata)
    end
end # function write_metadata

function prepare(
    step::Step{1},
    inputs::AbstractVector{<:AbstractString},
    template::PWscfInput,
    trial_eos::EquationOfState,
    pressures::AbstractVector{<:Real},
    metadatafiles::AbstractVector{<:AbstractString}
)
    # Checking parameters
    @assert length(inputs) == length(pressures) == length(metadatafiles) "The inputs, pressures and the metadata files must be the same size!"
    isnothing(template.cell_parameters) && (template = autofill_cell_parameters(template))
    if template.control.calculation != "scf"
        @warn "The calculation type is $(template.control.calculation), not \"scf\"! We will set it for you."
        @set! template.control.calculation = "scf"
    end
    # Write input and metadata files
    for (input, pressure) in zip(inputs, pressures)
        write_input(input, template, eos, pressure, verbose)
    end
    for (metadata, input) in zip(metadatafiles, inputs)
        write_metadata(metadata, template, input)
    end
end # function prepare
function prepare(
    step::Step{2},
    new_inputs::AbstractVector{<:AbstractString},
    previous_outputs::AbstractVector{<:AbstractString},
    template::PWscfInput,
    trial_eos::EquationOfState,
    pressures::AbstractVector{<:Real},
    metadatafiles::AbstractVector{<:AbstractString}
)
    # Checking parameters
    @assert length(new_inputs) == length(previous_outputs) == length(pressures) == length(metadatafiles) "The inputs, outputs, pressures and the metadata files must be the same size!"
    isnothing(template.cell_parameters) && (template = autofill_cell_parameters(template))
    if template.control.calculation != "vc-relax"
        @warn "The calculation type is $(template.control.calculation), not \"vc-relax\"! We will set it for you."
        @set! template.control.calculation = "vc-relax"
    end
    # Write input and metadata files
    for (input, pressure) in zip(new_inputs, pressures)
        write_input(input, template, trial_eos, pressure, verbose)
    end
    for (metadata, input) in zip(metadatafiles, inputs)
        write_metadata(metadata, template, input)
    end
end # function prepare

function finish(outputs::AbstractVector{<:AbstractString}, trial_eos::EquationOfState)
    energies = parse_total_energy.(outputs)
    volumes = prase_volume.(outputs)
    return lsqfit(EnergyForm(), trial_eos, volumes, energies)
end # function finish

# This is what the web service does
# function workflow(args)
#     prepare(Step(1), inputs, template, trial_eos, pressures, metadatafiles)
#     trial_eos = finish(outputs, trial_eos)
#     prepare(Step(2), inputs, previous_outputs, template, trial_eos, pressures, metadatafiles)
#     return finish(outputs, trial_eos)
# end # function workflow

end
