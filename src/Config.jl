module Config

using Configurations: @option

@option struct Directories
    root::String = pwd()
    prefix::String = "p="
    group_by_step::Bool = false
end

macro unit_vec_opt(type, unit, alias, criteria = (values, unit) -> nothing)
    return quote
        @option $alias struct $type
            values::AbstractVector
            unit::String
            function $type(values, unit = $unit)
                $criteria(values, unit)
                return new(values, unit)
            end
        end
        function $type(values::AbstractString, unit = $unit)
            values = eval(Meta.parse(values))
            typeassert(values, AbstractVector)
            return $type(values, unit)
        end
    end |> esc
end

end
