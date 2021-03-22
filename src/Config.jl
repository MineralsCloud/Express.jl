module Config

using Configurations: @option

@option struct Directories
    root::String = pwd()
    prefix::String = "p="
    group_by_step::Bool = false
end

macro vecunit(type, unit, alias = "", criteria = nothing)
    return quote
        @option $alias struct $type
            values::Any
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
