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
        error("unsupported eos name `\"$name\"`!")
    end
    parameters = config["parameters"]
    eos = if parameters isa AbstractVector
        ctor(map(myuparse, parameters))
    elseif parameters isa AbstractDict
        ctor(myuparse(parameters[string(f)]) for f in fieldnames(ctor))
    else
        error("unknown container of parameters!")
    end
    if !(eltype(eos) <: AbstractQuantity)
        @warn "the equation of state's elements seem not to have units! Be careful!"
    end
    return eos
end

function materialize_press(config)
    unit = myuparse(if haskey(config, "unit")
        config["unit"]
    else
        @info "no unit provided for `\"pressures\"`! \"GPa\" is assumed!"
        "GPa"
    end)
    return map(Base.Fix2(*, unit), config["values"])
end

function materialize_vol(config, templates)  # Arg `templates` have the same length as `pressures` already.
    if haskey(config, "volumes")
        subconfig = config["volumes"]
        unit = myuparse(if haskey(subconfig, "unit")
            subconfig["unit"]
        else
            @info "no unit provided for `\"volumes\"`! \"bohr^3\" is assumed!"
            "bohr^3"
        end)
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
    for key in ("np", "pressures", "bin", "templates")
        @assert haskey(config, key) "`\"$key\"` was not found in config!"
    end
    @assert config["np"] isa Integer && config["np"] >= 1
    checkconfig(currentsoftware(), config["bin"])  # To be implemented
    if haskey(config, "use_shell") && haskey(config, "shell_args") && config["use_shell"]
        @assert config["shell_args"] isa AbstractDict
    end
    let subconfig = config["pressures"], values = subconfig["values"]
        _alert(values)
        if length(config["templates"]) != 1
            if length(values) != length(config["templates"])
                throw(DimensionMismatch("templates and pressures have different lengths!"))
            end
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
