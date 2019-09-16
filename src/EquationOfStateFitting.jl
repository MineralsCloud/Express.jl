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
using QuantumESPRESSOParsers.OutputParsers.PWscf: read_total_energy,
                                                  read_cell_parameters,
                                                  read_head,
                                                  isjobdone
using Setfield: set, @lens

import ..Step
using ..SelfConsistentField: write_metadata

export update_alat_press, prepare, finish

function update_alat_press(template::PWscfInput, eos::EquationOfState, pressure::Real)
    volume = findvolume(PressureForm(), eos, pressure, (eps(), eos.v0 * 1.3))
    alat = cbrt(volume / det(template.cell_parameters.data))
    lenses = @batchlens(begin
        _.system.celldm ∘ _[$1]  # Get the `template`'s `system.celldm[1]` value
        _.cell.press             # Get the `template`'s `cell.press` value
    end)
    # Set the `template`'s `system.celldm[1]` and `cell.press` values with `alat` and `pressure`
    return set(template, lenses, (alat, pressure))
end # function update_alat_press

# This is a helper function and should not be exported.
function _validate(step::Step{N}, template::PWscfInput) where {N}
    control = @lens _.control
    lenses = @batchlens(
        control ∘ @lens _.calculation
        control ∘ @lens _.verbosity
        control ∘ @lens _.tstress
        control ∘ @lens _.tprnfor
    )
    template = set(template, lenses, (N == 1 ? "scf" : "vc-relax", "high", true, true))
    return isnothing(template.cell_parameters) ? autofill_cell_parameters(template) : template
end # function _validate

function prepare(
    step::Step,
    inputs::AbstractVector{<:AbstractString},
    template::PWscfInput,
    trial_eos::EquationOfState,
    pressures::AbstractVector{<:Real},
    metadatafiles::AbstractVector{<:AbstractString} = map(
        x -> splitext(x)[1] * ".json",
        inputs
    ),
    verbose::Bool = false
)
    # Check parameters
    @assert(
        length(inputs) == length(pressures) == length(metadatafiles),
        "The inputs, pressures and the metadata files must be the same size!"
    )
    template = _validate(step, template)
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
    ::Step{1},
    outputs::AbstractVector{<:AbstractString},
    trial_eos::EquationOfState
)
    energies = Float64[]
    volumes = Float64[]
    for output in outputs
        open(output, "r") do io
            s = read(io, String)
            push!(energies, (last ∘ read_total_energy)(s))
            push!(volumes, read_head(s)["unit-cell volume"])
        end
    end
    return lsqfit(EnergyForm(), trial_eos, volumes, energies)
end # function finish
function finish(
    ::Step{2},
    outputs::AbstractVector{<:AbstractString},
    trial_eos::EquationOfState
)
    energies = Float64[]
    volumes = Float64[]
    for output in outputs
        open(output, "r") do io
            s = read(io, String)
            isjobdone(s) || @warn("Job is not finished!")
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
