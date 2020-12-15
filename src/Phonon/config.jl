function materialize_press_vol(config)
    unit = uparse(
        if haskey(config, "unit")
            config["unit"]
        else
            @info "no unit provided for `\"pressures\"`! \"GPa\" is assumed!"
            u"GPa"
        end;
        unit_context = UNIT_CONTEXT,
    )
    return map(Base.Fix2(*, unit), config["values"])
end

function checkconfig(config)
    map(("templates", "qe", "workdir")) do key
        @assert haskey(config, key)
    end
    checkconfig(currentsoftware(), config["qe"])  # To be implemented
    let workdir = expanduser(config["workdir"])
        if !isdir(workdir)
            @warn "`\"workdir\"` \"$workdir\" is not reachable, be careful!"
        end
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
    subconfig = config[key]
    return
end

function materialize end
