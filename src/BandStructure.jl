"""
# module BandStructure



# Examples

```jldoctest
julia>
```
"""
module BandStructure

using Distances: euclidean
using QuantumESPRESSOBase.Cards.PWscf
using QuantumESPRESSOBase.Namelists.PWscf
using QuantumESPRESSOBase.Inputs.PWscf
using Setfield
using ShiftedArrays: circshift, lead

using Express
using Express.SelfConsistentField: write_metadata

export generate_path, update_kpoints, prepare

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
function generate_path(
    nodes::AbstractVector{<:AbstractVector},
    densities::AbstractVector{<:Integer} = 100 * ones(Int, length(nodes))
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

function update_kpoints(
    template::PWscfInput,
    nodes::AbstractVector{<:AbstractVector},
    densities::AbstractVector{<:Integer} = 100 * ones(Int, length(nodes))
)
    update_kpoints(template, generate_path(nodes, densities))
end # function update_kpoints
function update_kpoints(template::PWscfInput, path::AbstractVector{<:AbstractVector})
    data = map(x -> SpecialKPoint(x, 1), path)
    @set template.k_points = KPointsCard("crystal_b", data)
end # function update_kpoints

# This is a helper function and should not be exported
which_calculation(step::Step{1}) = "nscf"
which_calculation(step::Step{2}) = "bands"

# This is a helper function and should not be exported
function set_calculation(step::Step, template::PWscfInput)
    type = which_calculation(step)
    if template.control.calculation != calculation
        @warn "The calculation type is $(template.control.calculation), not \"$type\"! We will set it for you."
    end
    @set template.control.calculation = type  # Return a new `template` with its `control.calculation` to be `type`
end # function set_calculation

function prepare(
    step::Step{1},
    inputs::AbstractVector{<:AbstractString},
    template::PWscfInput,
    metadatafiles::AbstractVector{<:AbstractString}
)
    # Checking parameters
    @assert length(inputs) == length(metadatafiles) "The inputs and the metadata files must be the same size!"
    template = set_calculation(step, template)
    # Write input and metadata files
    for (input, metadata) in zip(inputs, metadatafiles)
        write_metadata(metadata, template, input)
    end
end # function prepare
function prepare(
    step::Step{2},
    inputs::AbstractVector{<:AbstractString},
    template::PWscfInput,
    nodes::AbstractVector{<:AbstractVector},
    densities::AbstractVector{<:Integer} = 100 * ones(Int, length(nodes))
    metadatafiles::AbstractVector{<:AbstractString}
)
    # Checking parameters
    @assert length(inputs) == length(metadatafiles) "The inputs and the metadata files must be the same size!"
    template = set_calculation(step, template)
    template = update_kpoints(template, nodes, densities)
    # Write input and metadata files
    for (input, metadata) in zip(inputs, metadatafiles)
        open(input, "r+") do io
            write(io, to_qe(object))
        end
        write_metadata(metadata, template, input)
    end
end # function prepare
function prepare(
    step::Step{3},
    inputs::AbstractVector{<:AbstractString},
    template::BandsNamelist,
    metadatafiles::AbstractVector{<:AbstractString}
)
    for (input, metadata) in zip(inputs, metadatafiles)
        open(input, "r+") do io
            write(io, to_qe(template))
        end
        write_metadata(metadata, template, input)
    end
end # function prepare

end
