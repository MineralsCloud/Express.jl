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
