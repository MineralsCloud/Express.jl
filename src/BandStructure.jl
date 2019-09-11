"""
# module BandStructure



# Examples

```jldoctest
julia>
```
"""
module BandStructure

using Distances: euclidean
using ShiftedArrays: circshift, lead

using Express
using Express.SelfConsistentField: write_metadata

export generate_path

abstract type PathType end
struct CircularPath <: PathType end
struct NoncircularPath <: PathType end

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
function generate_path(nodes::AbstractVector{<:AbstractVector}, densities::AbstractVector{<:Integer} = 100 * ones(Int, length(nodes)))
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

end
