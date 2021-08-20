module Config

using AbInitioSoftwareBase: load
using Configurations: @option

using ..Express: whichmodule

export loadconfig

@option struct Directories
    root::String = pwd()
    naming_convention::String = "p=%.1f"
    group_by_step::Bool = false
    Directories(root, naming_convention, group_by_step) =
        new(abspath(expanduser(root)), naming_convention, group_by_step)
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

function loadconfig(file)
    config = load(file)
    mod = whichmodule(pop!(config, "workflow"))
    return mod.Config.materialize(config)
end

end
