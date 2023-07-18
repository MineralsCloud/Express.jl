module Config

using Configurations: @option
using Formatting: sprintf1

@option struct InputFile
    name::String = "%s.in"
end

@option struct OutputFile
    name::String = "%s.out"
end

@option struct Directory
    root::String = pwd()
    name::String = "%s"
    input::InputFile = InputFile()
    output::OutputFile = OutputFile()
end

function getfiles(dir::Directory, name, filename)
    path = joinpath(dir.root, sprintf1(dir.name, name))
    input, output = sprintf1(dir.input.name, filename), sprintf1(dir.output.name, filename)
    return joinpath(path, input) => joinpath(path, output)
end

using Configurations: OptionField, option_m
using Unitful: Unitful, FreeUnits, Quantity, uparse, dimension, lookup_units
using UnitfulAtomic: UnitfulAtomic

import Configurations: from_dict

abstract type SamplingPoints end

function _uparse(str::AbstractString)
    return lookup_units([Unitful, UnitfulAtomic], Meta.parse(filter(!isspace, str)))
end

function from_dict(
    ::Type{<:SamplingPoints},
    ::OptionField{:numbers},
    ::Type{Vector{Float64}},
    str::AbstractString,
)
    return eval(Meta.parse(str))
end
function from_dict(
    ::Type{<:SamplingPoints}, ::OptionField{:unit}, ::Type{<:FreeUnits}, str::AbstractString
)
    return _uparse(str)
end

# Similar to https://github.com/JuliaCollections/IterTools.jl/blob/0ecaa88/src/IterTools.jl#L1028-L1032
function Base.iterate(iter::SamplingPoints, state=1)
    if state > length(iter.numbers)
        return nothing
    else
        return getindex(iter.numbers, state) * iter.unit, state + 1
    end
end

Base.eltype(iter::SamplingPoints) = Quantity{Float64,dimension(iter.unit),typeof(iter.unit)}

Base.length(iter::SamplingPoints) = length(iter.numbers)

Base.size(iter::SamplingPoints) = size(iter.numbers)

end
