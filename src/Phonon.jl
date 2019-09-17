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
function _preset(template::PhononInput)
    # TODO: Implement this
end # function _preset

# This is a helper function and should not be exported.
"""
    _inject_shared_info(a::PhononInput, b::PWscfInput)

Inject shared information from a `PWscfInput` to a `PhononInput`.

A `PWscfInput` before a `PhononInput` has the information of `outdir` and `prefix`. They must keep the same in a
phonon calculation.
"""
function _inject_shared_info(a::PhononInput, b::PWscfInput)
    # TODO: Implement this
    control = @lens _.control
    lenses = @batchlens(begin
        control ∘ @lens _.outdir
        control ∘ @lens _.prefix
    end)
    newlenses = @batchlens(begin
        _.phonon.outdir
        _.phonon.prefix
    end)
    return set(a, newlenses, get(b, lenses))
end # function _inject_shared_info

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
    # Check parameters
    @assert(
        length(inputs) == length(outputs) == length(metadatafiles),
        "The inputs, outputs and the metadata files must be the same length!"
    )
    template = _preset(template)
    for (input, output, structure, metadatafile) in zip(inputs, outputs, metadatafiles)
        # Get a new `object` from the `template`, with its `alat` and `pressure` changed
        object = update_structure(output, template)
        write(input, to_qe(object, verbose = verbose))  # Write the `object` to a Quantum ESPRESSO input file
        write_metadata(metadatafile, input, template)
    end
end # function prepare
function prepare(
    step::Step{2},
    phonon_inputs::AbstractVector{<:AbstractString},
    pwscf_inputs::AbstractVector{<:AbstractString},
    template::PhononInput,
    verbose::Bool = false
)
    # Check parameters
    @assert(
        length(phonon_inputs) == length(pwscf_inputs),
        "The PWscf and the PHonon inputs files must have the same length!"
    )
    template = _preset(template)
    for (phonon_input, pwscf_input) in zip(phonon_inputs, pwscf_inputs)
        open(pwscf_input, "r") do io
            object = parse(PWscfInput, read(io, String))
        end
        template = _inject_shared_info(template, object)
        write(phonon_input, template)
    end
end # function prepare

end
