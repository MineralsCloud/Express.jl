function materialize_eos(config)
    name = lowercase(config["name"])
    ctor = if name in ("m", "murnaghan")
        Murnaghan
    elseif name in ("bm2", "birchmurnaghan2nd", "birch-murnaghan-2")
        BirchMurnaghan2nd
    elseif name in ("bm3", "birchmurnaghan3rd", "birch-murnaghan-3")
        BirchMurnaghan3rd
    elseif name in ("bm4", "birchmurnaghan4th", "birch-murnaghan-4")
        BirchMurnaghan4th
    elseif name in ("pt2", "poiriertarantola2nd", "poirier-tarantola-2")
        PoirierTarantola2nd
    elseif name in ("pt3", "poiriertarantola3rd", "poirier-tarantola-3")
        PoirierTarantola3rd
    elseif name in ("pt4", "poiriertarantola4th", "poirier-tarantola-4")
        PoirierTarantola4th
    elseif name in ("v", "vinet")
        Vinet
    else
        error("unsupported eos type `\"$type\"`!")
    end
    values = (
        v * uparse(string(u); unit_context = UNIT_CONTEXT) for
        (v, u) in config["parameters"]
    )
    return ctor(values...)
end

function materialize_press(config)
    unit = uparse(
        if haskey(config, "unit")
            config["unit"]
        else
            @info "no unit provided for `\"pressures\"`! \"GPa\" is assumed!"
            "GPa"
        end;
        unit_context = UNIT_CONTEXT,
    )
    return map(Base.Fix2(*, unit), config["values"])
end

function materialize_vol(config, templates)  # Arg `templates` have the same length as `pressures` already.
    if haskey(config, "volumes")
        subconfig = config["volumes"]
        unit = uparse(
            if haskey(subconfig, "unit")
                subconfig["unit"]
            else
                @info "no unit provided for `\"volumes\"`! \"bohr^3\" is assumed!"
                "bohr^3"
            end;
            unit_context = UNIT_CONTEXT,
        )
        if length(subconfig["values"]) == 1
            return repeat(subconfig["values"] * unit, length(templates))
        else
            return map(Base.Fix2(*, unit), subconfig["values"])
        end
    else
        return map(cellvolume, templates) * bohr^3  # FIXME: Are units all `bohr^3` for different software?
    end
end

function materialize_dirs(config, pressures)
    return map(pressures) do pressure
        abspath(joinpath(expanduser(config), "p" * string(ustrip(pressure))))
    end
end

function _alert(pressures)
    if length(pressures) <= 5
        @info "pressures <= 5 may give unreliable results, consider more if possible!"
    end
    if minimum(pressures) >= zero(eltype(pressures))
        @warn "for better fitting, we need at least 1 negative pressure!"
    end
end

function checkconfig(config)
    for key in ("pressures", "qe", "templates", "workdir")
        @assert haskey(config, key) "`\"$key\"` was not found in config!"
    end
    checkconfig(currentsoftware(), config["qe"])  # To be implemented
    let subconfig = config["pressures"], values = subconfig["values"]
        _alert(values)
        if length(config["templates"]) != 1
            if length(values) != length(config["templates"])
                throw(DimensionMismatch("templates and pressures have different lengths!"))
            end
        end
    end
    let workdir = expanduser(config["workdir"])
        if !isdir(workdir)
            @warn "`\"workdir\"` \"$workdir\" is not reachable, be careful!"
        end
    end
    for path in config["templates"]
        if !isfile(path)
            @warn "template \"$path\" is not reachable, be careful!"
        end
    end
    if haskey(config, "trial_eos")
        @assert !haskey(config, "volumes") "key \"trial_eos\" and \"volumes\" are mutually exclusive!"
        for key in ("name", "parameters")
            @assert haskey(config["trial_eos"], key) "the trial eos needs `\"$key\"` specified!"
        end
    end
    if haskey(config, "volumes")
        subconfig = config["volumes"]
        if length(subconfig["values"]) != 1
            if length(subconfig["values"]) != length(config["pressures"]["values"])
                throw(DimensionMismatch("volumes and pressures have different lengths!"))
            end
        end
    end
    return
end

function materialize end
