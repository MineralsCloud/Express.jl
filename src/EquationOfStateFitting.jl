"""
# module EquationOfStateFitting



# Examples

```jldoctest
julia>
```
"""
module EquationOfStateFitting

using LinearAlgebra: det

using Compat: isnothing
using EquationsOfState
using EquationsOfState.Collections: EquationOfState
using EquationsOfState.NonlinearFitting: lsqfit
using EquationsOfState.FindVolume: findvolume
using Kaleido: @batchlens
using QuantumESPRESSOBase: to_qe
using QuantumESPRESSOBase.Inputs.PWscf: PWscfInput
using QuantumESPRESSOParsers.OutputParsers.PWscf: read_total_energy, read_cell_parameters
using Setfield: set, @lens

using Express: Step
using Express.SelfConsistentField: write_metadata

export update_alat_press, write_input, prepare, finish

function update_alat_press(template::PWscfInput, eos::EquationOfState, pressure::Real)
    volume = findvolume(PressureForm(), eos, pressure, (0, 1000), Order8())
    alat = cbrt(volume / det(template.cell_parameters.data))
    lenses = @batchlens(begin
        _.system.celldm ∘ _[$1]  # Get the `template`'s `system.celldm[1]` value
        _.cell.press             # Get the `template`'s `cell.press` value
    end)
    # Set the `template`'s `system.celldm[1]` and `cell.press` values with `alat` and `pressure`
    return set(template, lenses, (alat, pressure))
end # function update_alat_press

function write_input(
    input::AbstractString,
    template::PWscfInput,
    eos::EquationOfState,
    pressure::Real,
    verbose::Bool = false
)
    # Get a new `object` from the `template`, with its `alat` and `pressure` changed
    object = update_alat_press(template, eos, pressure)
    # Write the `object` to a Quantum ESPRESSO input file
    open(input, "r+") do io
        write(io, to_qe(object, verbose = verbose))
    end
    return
end # function write_input

# This is a helper function and should not be exported
_calculationof(step::Step{1}) = "scf"
_calculationof(step::Step{2}) = "vc-relax"

# This is a helper function and should not be exported
function _set_calculation(step::Step, template::PWscfInput)
    type = _calculationof(step)
    lens = @lens _.control.calculation
    if get(template, lens) != type
        @warn("The calculation type of step $step should be \"$type\", not $(get(template, lens))! I will set it for you.")
    end
    return set(template, lens, type)  # Return a new `template` whose `control.calculation` is `type`
end # function _set_calculation

function prepare(
    step::Step,
    inputs::AbstractVector{<:AbstractString},
    template::PWscfInput,
    trial_eos::EquationOfState,
    pressures::AbstractVector{<:Real},
    metadatafiles::AbstractVector{<:AbstractString} = map(x -> splitext(x)[1] * ".json", inputs),
    verbose::Bool = false
)
    # Checking parameters
    @assert length(inputs) == length(pressures) == length(metadatafiles) "The inputs, pressures and the metadata files must be the same size!"
    isnothing(template.cell_parameters) && (template = autofill_cell_parameters(template))
    template = _set_calculation(step, template)
    # Write input and metadata files
    for (input, pressure, metadata) in zip(inputs, pressures, metadatafiles)
        write_input(input, template, trial_eos, pressure, verbose)
        write_metadata(metadata, template, input)
    end
    return
end # function prepare

function finish(outputs::AbstractVector{<:AbstractString}, trial_eos::EquationOfState)
    energies = Float64[]
    volumes = Float64[]
    for output in outputs
        open(output, "r") do io
            s = read(io, String)
            isjobdone(s) || @warn "Job is not finished!"
            push!(energies, (last ∘ read_total_energy)(s))
            push!(volumes, (det ∘ last ∘ read_cell_parameters)(s))
        end
    end
    return lsqfit(EnergyForm(), trial_eos, volumes, energies)
end # function finish

# This is what the web service does
# function workflow(args)
#     prepare(Step(1), inputs, template, trial_eos, pressures, metadatafiles)
#     trial_eos = finish(outputs, trial_eos)
#     prepare(Step(2), inputs, template, trial_eos, pressures, metadatafiles)
#     return finish(outputs, trial_eos)
# end # function workflow

end
