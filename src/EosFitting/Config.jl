module Config

using Compat: isnothing
using Configurations: from_dict, @option
using Crystallography: cellvolume
using EquationsOfStateOfSolids:
    Murnaghan,
    BirchMurnaghan2nd,
    BirchMurnaghan3rd,
    BirchMurnaghan4th,
    PoirierTarantola2nd,
    PoirierTarantola3rd,
    PoirierTarantola4th,
    Vinet
using Unitful: AbstractQuantity, ustrip

using ...Express: myuparse

export materialize_eos

@option struct Templates
    paths::AbstractVector
    function Templates(paths)
        @assert length(paths) >= 1
        for path in paths
            if !isfile(path)
                @warn "template \"$path\" is not reachable, be careful!"
            end
        end
        return new(paths)
    end
end

@option "pressures" struct Pressures
    values::AbstractVector
    unit::String = "GPa"
    function Pressures(values, unit)
        if length(values) <= 5
            @info "less than 6 pressures may not fit accurately, consider adding more!"
        end
        if minimum(values) >= zero(eltype(values))
            @warn "for better fitting result, provide at least 1 negative pressure!"
        end
        return new(values, unit)
    end
end

@option "volumes" struct Volumes
    values::AbstractVector
    unit::String = "bohr^3"
end

@option "trial_eos" struct TrialEos
    name::String
    parameters::AbstractVector
end

@option "fit" struct EosFittingConfig
    templates::Templates
    fixed::Union{Pressures,Volumes,Nothing}
    trial_eos::Union{TrialEos,Nothing}
    function EosFittingConfig(templates, fixed, trial_eos)
        if length(templates.paths) != 1  # Always >= 1
            if length(templates.paths) != length(fixed.values)
                throw(
                    DimensionMismatch(
                        "templates and pressures or volumes have different lengths!",
                    ),
                )
            end
        end
        return new(templates, fixed, trial_eos)
    end
end

function materialize_eos(config::TrialEos)
    name = config.name
    T = if name in ("m", "murnaghan")
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
    eos = _materialize_eos(T, config.parameters)
    if !(eltype(eos) <: AbstractQuantity)
        @warn "the equation of state's elements seem not to have units! Be careful!"
    end
    return eos
end
materialize_eos(config::AbstractDict) = materialize_eos(from_dict(TrialEos, config))
_materialize_eos(T, parameters::AbstractVector) = T(map(myuparse, parameters))
_materialize_eos(T, parameters::AbstractDict) =
    T((myuparse(parameters[string(f)]) for f in fieldnames(T))...)

function materialize_press(config::Pressures)
    unit = myuparse(config.unit)
    return config.values .* unit
end
materialize_press(config::AbstractDict) = materialize_press(from_dict(Pressures, config))

function materialize_vol(config::Volumes)
    unit = myuparse(config.unit)
    return config.values .* config.unit
end
function materialize_vol(config::EosFittingConfig)
    if isnothing(config.fixed)  # If no volume or pressure is provided, use templates cell volumes
        return materialize_vol(Volumes(map(cellvolume, config.templates)))
    elseif config.fixed isa Pressures
        error("wrong method called! Use `materialize_press` instead!")
    else
        return materialize_vol(config.fixed)
    end
end
materialize_vol(config::AbstractDict) = materialize_vol(from_dict(EosFittingConfig, config))

function materialize_dirs(config, pressures)
    return map(pressures) do pressure
        abspath(joinpath(expanduser(config), "p" * string(ustrip(pressure))))
    end
end

function checkconfig(config)
    for key in ("np", "pressures", "bin", "templates")
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
