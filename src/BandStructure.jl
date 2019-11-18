"""
# module BandStructure



# Examples

```jldoctest
julia>
```
"""
module BandStructure

using QuantumESPRESSO: to_qe
using QuantumESPRESSO.Cards.PWscf: SpecialKPoint, KPointsCard
using QuantumESPRESSO.Namelists.PWscf: BandsNamelist
using QuantumESPRESSO.Inputs.PWscf: PWInput
using Setfield: @set
using ShiftedArrays: circshift, lead

import ..Step
using ..SelfConsistentField: write_metadata

export generate_path, update_kpoints, prepare

abstract type PathStyle end
struct CircularPath <: PathStyle end
struct NoncircularPath <: PathStyle end

# This is a helper function and should not be exported
euclidean(x, y) = sqrt(sum((x - y) .^ 2))

"""
    generate_path(nodes, densities = 100 * ones(Int, length(nodes)))

Generate a reciprocal space path from each node.

# Arguments
- `nodes::AbstractVector{<:AbstractVector}`: a vector of 3-element vectors.
- `densities::AbstractVector{<:Integer}`: the default value is a circular path.

# Examples
```jldoctest
julia> nodes = [
    [0.0, 0.0, 0.5],
    [0.0, 0.0, 0.0],
    [0.3333333333, 0.3333333333, 0.5],
    [0.3333333333, 0.3333333333, -0.5],
    [0.3333333333, 0.3333333333, 0.0],
    [0.5, 0.0, 0.5],
    [0.5, 0.0, 0.0]
];

julia> BandStructure.generate_path(nodes)  # Generate a circular path
693-element Array{Any,1}:
...

julia> BandStructure.generate_path(nodes, 100 * ones(Int, length(nodes) - 1))  # Generate a noncircular path
594-element Array{Any,1}:
...
```
"""
function generate_path(
    nodes::AbstractVector{<:AbstractVector},
    densities::AbstractVector{<:Integer} = 100 * ones(Int, length(nodes)),
)
    if length(densities) == length(nodes)
        _generate_path(nodes, densities, CircularPath())
    elseif length(densities) == length(nodes) - 1
        _generate_path(nodes, densities, NoncircularPath())
    else
        error("The length of `densities` should be either length of `length(nodes)` or `length(nodes) - 1`!")
    end
end # function generate_path
function _generate_path(nodes, densities, ::CircularPath)
    path = []
    for (thisnode, nextnode, density) in zip(nodes, circshift(nodes, -1), densities)
        distance = euclidean(thisnode, nextnode)  # Compute Euclidean distance between two vectors
        step = (nextnode - thisnode) / distance
        for x in range(0, stop = distance, length = density - 1)
            push!(path, thisnode + x * step)
        end
    end
    return path
end # function _generate_path
function _generate_path(nodes, densities, ::NoncircularPath)
    path = []
    for (thisnode, nextnode, density) in zip(nodes, lead(nodes), densities)
        ismissing(nextnode) && break
        distance = euclidean(thisnode, nextnode)  # Compute Euclidean distance between two vectors
        step = (nextnode - thisnode) / distance
        for x in range(0, stop = distance, length = density - 1)
            push!(path, thisnode + x * step)
        end
    end
    return path
end # function _generate_path

function update_kpoints(template::PWInput, path::AbstractVector{<:AbstractVector})
    data = map(x -> SpecialKPoint(x, 1), path)
    @set template.k_points = KPointsCard("crystal_b", data)
end # function update_kpoints

# This is a helper function and should not be exported
which_calculation(step::Step{1}) = "nscf"
which_calculation(step::Step{2}) = "bands"

# This is a helper function and should not be exported
function set_calculation(step::Step, template::PWInput)
    type = which_calculation(step)
    if template.control.calculation != type
        @warn "The calculation type is $(template.control.calculation), not \"$type\"! We will set it for you."
    end
    @set template.control.calculation = type  # Return a new `template` with its `control.calculation` to be `type`
end # function set_calculation

function prepare(
    step::Step{1},
    inputs::AbstractVector{<:AbstractString},
    template::PWInput,
    metadatafiles::AbstractVector{<:AbstractString},
)
    # Checking parameters
    @assert length(inputs) == length(metadatafiles) "The inputs and the metadata files must be the same size!"
    template = set_calculation(step, template)
    # Write input and metadata files
    for (input, metadata) in zip(inputs, metadatafiles)
        write_metadata(metadata, input, template)
    end
end # function prepare
function prepare(
    step::Step{2},
    inputs::AbstractVector{<:AbstractString},
    template::PWInput,
    nodes::AbstractVector{<:AbstractVector},
    densities::AbstractVector{<:Integer} = 100 * ones(Int, length(nodes)),
    metadatafiles::AbstractVector{<:AbstractString} = map(
        x -> splitext(x)[1] * ".json",
        inputs
    )
)
    # Checking parameters
    @assert length(inputs) == length(metadatafiles) "The inputs and the metadata files must be the same size!"
    template = set_calculation(step, template)
    template = update_kpoints(template, generate_path(nodes, densities))
    # Write input and metadata files
    for (input, metadata) in zip(inputs, metadatafiles)
        open(input, "r+") do io
            write(io, to_qe(template))
        end
        write_metadata(metadata, input, template)
    end
end# function prepare
function prepare(
    step::Step{3},
    inputs::AbstractVector{<:AbstractString},
    template::BandsNamelist,
    metadatafiles::AbstractVector{<:AbstractString},
)
    for (input, metadata) in zip(inputs, metadatafiles)
        open(input, "r+") do io
            write(io, to_qe(template))
        end
        write_metadata(metadata, input, template)
    end
end # function prepare

end
