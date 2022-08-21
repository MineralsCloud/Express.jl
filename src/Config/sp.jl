using Configurations: OptionField, option_m
using Unitful: Unitful, FreeUnits, Quantity, uparse, dimension, lookup_units
using UnitfulAtomic: UnitfulAtomic

import Configurations: from_dict

abstract type SamplingPoints end

# See https://github.com/Roger-luo/Configurations.jl/blob/933fd46/src/codegen.jl#L82-L84
macro sp(type, unit, alias, check=(_, _) -> nothing)
    unit = _uparse(unit)
    ex = :(struct $type <: $SamplingPoints
        numbers::Vector{Float64}
        unit::$FreeUnits
        function $type(numbers, unit=$unit)
            $check(numbers, unit)
            return new(numbers, unit)
        end
    end)
    return esc(option_m(__module__, ex, alias))
end

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
