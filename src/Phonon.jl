"""
# module Phonon



# Examples

```jldoctest
julia>
```
"""
module Phonon

using Kaleido: @batchlens
using QuantumESPRESSOBase: to_qe
using QuantumESPRESSOBase.Inputs.PWscf: PWscfInput
using QuantumESPRESSOParsers.OutputParsers.PWscf: read_cell_parameters, read_atomic_positions
using Setfield: set

import ..Step
using ..SelfConsistentField: write_metadata

export update_structure, prepare

"""
    update_structure(output::AbstractString, template::PWscfInput)

Read structure information from `output`, and update related fields of `template`.
"""
function update_structure(output::AbstractString, template::PWscfInput)
    cell_parameters = read_cell_parameters(output)
    atomic_positions = read_atomic_positions(output)  # TODO: Implement `read_atomic_positions`
    lenses = @batchlens(begin
        _.system.celldm ∘ _[$1]
        _.atomic_positions
        _.cell_parameters
    end)
    return set(template, lenses, (1, atomic_positions, cell_parameters))
end # function update_structure

# This is a helper function and should not be exported.
function _preset(template::PWscfInput)
    control = @lens _.control
    lenses = @batchlens(begin
        control ∘ @lens _.calculation  # Get the `template`'s `control.calculation` value
        control ∘ @lens _.verbosity    # Get the `template`'s `control.verbosity` value
        control ∘ @lens _.tstress      # Get the `template`'s `control.tstress` value
        control ∘ @lens _.tprnfor      # Get the `template`'s `control.tprnfor` value
    end)
    # Set the `template`'s values with...
    template = set(template, lenses, ("scf", "high", true, true))
    return isnothing(template.cell_parameters) ? autofill_cell_parameters(template) : template
end # function _preset

"""
    prepare(step::Step{1}, inputs, outputs, template, metadatafiles[, verbose::Bool = false])

Prepare input files of the first step of a phonon calculation.

# Arguments
- `step::Step{1}`: Denotes the first step of a phonon calculation.
- `inputs::AbstractVector{<:AbstractString}`: The input files of Quantum ESPRESSO of a phonon calculation. If not exist,
   they will be automatically created.
- `outputs::AbstractVector{<:AbstractString}`: The output files of Quantum ESPRESSO of a vc-relax calculation.
- `template::PWscfInput`:
- `metadatafiles::AbstractVector{<:AbstractString}`:
- `verbose::Bool = false`: control the format of input files, verbose or not.
"""
function prepare(
    step::Step{1},
    inputs::AbstractVector{<:AbstractString},
    outputs::AbstractVector{<:AbstractString},
    template::PWscfInput,
    metadatafiles::AbstractVector{<:AbstractString} = map(x -> splitext(x)[1] * ".json", inputs),
    verbose::Bool = false
)
    template = _preset(template)
    for (input, output, structure, metadatafile) in zip(inputs, outputs, metadatafiles)
        # Get a new `object` from the `template`, with its `alat` and `pressure` changed
        object = update_structure(output, template)
        write(input, to_qe(object, verbose = verbose))  # Write the `object` to a Quantum ESPRESSO input file
        write_metadata(metadatafile, input, template)
    end
end # function prepare

end
