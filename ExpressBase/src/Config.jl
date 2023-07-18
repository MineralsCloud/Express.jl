module Config

using Configurations: OptionField
using Unitful: Unitful, FreeUnits, Quantity, uparse, dimension, lookup_units
using UnitfulAtomic

import Configurations: from_dict

export Subdirectory, list_io_files

using Configurations: @option
using Formatting: sprintf1

@option struct InputFile
    base::String = ""
    extension::String = "in"
end

@option struct OutputFile
    base::String = ""
    extension::String = "out"
end

@option struct Subdirectory
    root::String = pwd()
    pattern::String = "%s"

@option "io" struct IO
    subdir::Subdirectory = Subdirectory()
    in::InputFile = InputFile()
    out::OutputFile = OutputFile()
end

function list_io_files(dir::Subdirectory, name, filename)
    path = joinpath(dir.root, sprintf1(dir.name, name))
    input, output = sprintf1(dir.input.name, filename), sprintf1(dir.output.name, filename)
    return joinpath(path, input) => joinpath(path, output)
end

abstract type SamplingPoints end

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

_uparse(str::AbstractString) =
    lookup_units([Unitful, UnitfulAtomic], Meta.parse(filter(!isspace, str)))

end
