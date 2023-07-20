module Config

using Configurations: OptionField, @option
using Formatting: sprintf1
using Unitful: Unitful, FreeUnits, Quantity, uparse, dimension, lookup_units
using UnitfulAtomic

import Configurations: from_dict

abstract type AbstractConfig end

@option struct InputFile <: AbstractConfig
    extension::String = "in"
end

@option struct OutputFile <: AbstractConfig
    extension::String = "out"
end

@option struct Subdirectory <: AbstractConfig
    root::String = pwd()
    pattern::String = "%s"
end

@option struct IO <: AbstractConfig
    subdir::Subdirectory = Subdirectory()
    in::InputFile = InputFile()
    out::OutputFile = OutputFile()
end

function list_io(io::IO, dir, file)
    path = joinpath(io.subdir.root, sprintf1(io.subdir.pattern, dir))
    in, out = join((file, io.in.extension), '.'), join((file, io.out.extension), '.')
    return joinpath(path, in) => joinpath(path, out)
end

abstract type SamplingPoints <: AbstractConfig end

from_dict(
    ::Type{<:SamplingPoints},
    ::OptionField{:numbers},
    ::Type{Vector{Float64}},
    str::AbstractString,
) = eval(Meta.parse(str))
from_dict(
    ::Type{<:SamplingPoints}, ::OptionField{:unit}, ::Type{<:FreeUnits}, str::AbstractString
) = _uparse(str)

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

_uparse(str::AbstractString) =
    lookup_units([Unitful, UnitfulAtomic], Meta.parse(filter(!isspace, str)))

end
