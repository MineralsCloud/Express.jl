using Configurations: OptionField, option_m
using Unitful: Unitful, FreeUnits, Quantity, uparse, dimension, lookup_units
using UnitfulAtomic: UnitfulAtomic

import Configurations: from_dict

export @vopt

abstract type VectorWithUnitOption end

# See https://github.com/Roger-luo/Configurations.jl/blob/933fd46/src/codegen.jl#L82-L84
macro vopt(type, unit, alias, check = (_, _) -> nothing)
    unit = _uparse(unit)
    ex = :(struct $type <: $VectorWithUnitOption
        vector::Vector{Float64}
        unit::$FreeUnits
        function $type(vector, unit = $unit)
            $check(vector, unit)
            return new(vector, unit)
        end
    end)
    return esc(option_m(__module__, ex, alias))
end

_uparse(str::AbstractString) =
    lookup_units([Unitful, UnitfulAtomic], Meta.parse(filter(!isspace, str)))

from_dict(
    ::Type{<:VectorWithUnitOption},
    ::OptionField{:vector},
    ::Type{Vector{Float64}},
    str::AbstractString,
) = eval(Meta.parse(str))
from_dict(
    ::Type{<:VectorWithUnitOption},
    ::OptionField{:unit},
    ::Type{<:FreeUnits},
    str::AbstractString,
) = _uparse(str)

# Similar to https://github.com/JuliaCollections/IterTools.jl/blob/0ecaa88/src/IterTools.jl#L1028-L1032
function Base.iterate(iter::VectorWithUnitOption, state = 1)
    if state > length(iter.vector)
        return nothing
    else
        return getindex(iter.vector, state) * iter.unit, state + 1
    end
end

Base.eltype(iter::VectorWithUnitOption) =
    Quantity{Float64,dimension(iter.unit),typeof(iter.unit)}

Base.length(iter::VectorWithUnitOption) = length(iter.vector)

Base.size(iter::VectorWithUnitOption) = size(iter.vector)
