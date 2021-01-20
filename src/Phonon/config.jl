function materialize_press_vol(config)
    unit = myuparse(if haskey(config, "unit")
        config["unit"]
    else
        @info "no unit provided for `\"pressures\"`! \"GPa\" is assumed!"
        u"GPa"
    end)
    return map(Base.Fix2(*, unit), config["values"])
end

function checkconfig(config)
    for key in ("np", "bin", "templates")
        @assert haskey(config, key) "`\"$key\"` was not found in config!"
    end
    @assert config["np"] isa Integer && config["np"] >= 1
    checkconfig(currentsoftware(), config["bin"])  # To be implemented
    if haskey(config, "use_shell") && config["use_shell"]
        @assert config["shell_args"] isa AbstractDict
    end
    for paths in config["templates"]
        for path in paths
            if !isfile(path)
                @warn "template \"$path\" is not reachable, be careful!"
            end
        end
    end
    if haskey(config, "pressures")
        key = "pressures"
    elseif haskey(config, "volumes")
        key = "volumes"
    else
        error("\"pressures\" or \"volumes\", there must be one!")
    end
    return
end

function materialize end
