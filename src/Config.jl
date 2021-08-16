module Config

using Configurations: @option

@option struct Directories
    root::String = pwd()
    naming_convention::String = "p=%.1f"
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
        $type(values::AbstractString, unit = $unit) =
            $type(eval(Meta.parse(values)), unit)
    end |> esc
end

end
