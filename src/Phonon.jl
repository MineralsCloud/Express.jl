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
using Setfield: set

import ..Step
using ..SelfConsistentField: write_metadata

export update_structure, prepare

function update_structure(output::AbstractString, template::PWscfInput)
    structure = readvolume(outputs)
    lenses = @batchlens begin
        _.system.celldm ∘ _[$1]
        _.cell_parameters
    end
    set(template, lenses, (structure[:celldm], structure[:cell_parameters]))
end # function update_structure

function prepare(
    step::Step{1},
    inputs::AbstractVector{<:AbstractString},
    template::PWscfInput,
    pressures::AbstractVector{<:Real},
    metadatafiles::AbstractVector{<:AbstractString}
)
    objects = update_structure.(outputs, template)
    for input in inputs
        open(input, "r+") do io
            write(io, to_qe(object))
        end
    end
    for (metadata, input) in zip(metadatafiles, inputs)
        write_metadata(metadata, input, template)
    end
end # function prepare

end
