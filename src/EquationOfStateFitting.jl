"""
# module EquationOfStateFitting



# Examples

```jldoctest
julia>
```
"""
module EquationOfStateFitting

using LinearAlgebra

using EquationsOfState
using EquationsOfState.Collections
using EquationsOfState.NonlinearFitting
using EquationsOfState.FindVolume
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
using Express.SelfConsistentField: write_metadata

export update_alat_press, write_input, prepare, finish

function update_alat_press(template::PWscfInput, eos::EquationOfState, pressure::Real)
    volume = findvolume(PressureForm(), eos, pressure, (0, 1000), Order8())
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
    @assert length(new_inputs) == length(previous_outputs) == length(pressures) ==
            length(metadatafiles) "The inputs, outputs, pressures and the metadata files must be the same size!"
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
