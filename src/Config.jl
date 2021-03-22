module Config

using Configurations: @option

macro vecunit(type, unit, alias = "", criteria = nothing)
    return quote
        @option $alias struct $type
            values
            unit::String
            function ($type)(values::AbstractVector, unit = $unit)
                $criteria
                return new(values, unit)
            end
        end
        function ($type)(values::AbstractString, unit = $unit)
            values = eval(Meta.parse(values))
            typeassert(values, AbstractVector)
            return ($type)(values, unit)
        end
    end |> esc
end

end
