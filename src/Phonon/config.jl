function checkconfig(settings)
    map(("templates", "pressures", "workdir")) do key
        @assert haskey(settings, key)
    end
    if !isdir(expanduser(settings["workdir"]))
        @warn "`workdir` is not reachable, be careful!"
    end
    for paths in settings["templates"]
        for path in paths
            if !isfile(path)
                @warn "template \"$path\" is not reachable, be careful!"
            end
        end
    end
end
