module Config

using Configurations: @option

using ...Express: myuparse

@option struct Templates
    paths::AbstractVector
    function Templates(paths)
        @assert length(paths) >= 1
        for path in paths
            if !isfile(path)
                @warn "template \"$path\" is not reachable, be careful!"
            end
        end
        return Templates(paths)
    end
end

@option "pressures" struct Pressures
    values::AbstractVector
    unit::String = "GPa"
end

@option "volumes" struct Volumes
    values::AbstractVector
    unit::String = "bohr^3"
end

@option "fit" struct PhononConfig
    templates::Templates
    fixed::Union{Pressures,Volumes}
    function PhononConfig(templates, fixed, trial_eos)
        if length(templates.paths) != 1  # Always >= 1
            if length(templates.paths) != length(fixed.values)
                throw(
                    DimensionMismatch(
                        "templates and pressures or volumes have different lengths!",
                    ),
                )
            end
        end
        return new(templates, fixed)
    end
end

function materialize_press_vol(config::Union{Pressures,Volumes})
    unit = myuparse(config.unit)
    return config.values .* unit
end

function checkconfig(config)
    for key in ("np", "bin", "templates")
        @assert haskey(config, key) "`\"$key\"` was not found in config!"
    end
    checkconfig(currentsoftware(), config["bin"])  # To be implemented
    if haskey(config, "use_shell") && haskey(config, "shell_args") && config["use_shell"]
        @assert config["shell_args"] isa AbstractDict
    end
    return
end

function materialize end

end
