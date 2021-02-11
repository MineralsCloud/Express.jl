module Config

using AbInitioSoftwareBase.Cli: CliConfig
using Compat: isnothing
using Configurations: @option
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
using EquationsOfStateOfSolids.Inverse: NumericalInversionOptions
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

@option struct TrialEos
    name::String
    parameters::Union{AbstractVector,AbstractDict}
end

@option struct Directories
    root::String = pwd()
    prefix::String = "p="
    group_by_step::Bool = false
end

@option struct EosFittingConfig{T<:CliConfig}
    templates::Templates
    trial_eos::TrialEos
    fixed::Union{Pressures,Volumes,Nothing} = nothing
    dirs::Directories = Directories()
    num_inv::NumericalInversionOptions = NumericalInversionOptions()
    cli::T
    function EosFittingConfig{T}(
        templates,
        trial_eos,
        fixed,
        dirs,
        num_inv,
        cli::T,
    ) where {T}
        if length(templates.paths) != 1  # Always >= 1
            if !isnothing(fixed)
                if length(templates.paths) != length(fixed.values)
                    throw(
                        DimensionMismatch(
                            "templates and pressures or volumes have different lengths!",
                        ),
                    )
                end
            end
        end
        return new(templates, trial_eos, fixed, dirs, num_inv, cli)
    end
end

function materialize_eos(config::TrialEos)
    name = filter(c -> isletter(c) || isdigit(c), lowercase(config.name))
    T = if name == "m" || occursin("murnaghan", name)
        Murnaghan
    elseif name == "bm2" || occursin("birchmurnaghan2", name)
        BirchMurnaghan2nd
    elseif name == "bm3" || occursin("birchmurnaghan3", name)
        BirchMurnaghan3rd
    elseif name == "bm4" || occursin("birchmurnaghan4", name)
        BirchMurnaghan4th
    elseif name == "pt2" || occursin("poiriertarantola2", name)
        PoirierTarantola2nd
    elseif name == "pt3" || occursin("poiriertarantola3", name)
        PoirierTarantola3rd
    elseif name == "pt4" || occursin("poiriertarantola4", name)
        PoirierTarantola4th
    elseif name == "v" || occursin("vinet", name)
        Vinet
    else
        error("unsupported eos name `\"$name\"`!")
    end
    return _materialize_eos(T, config.parameters)
end
_materialize_eos(T, parameters::AbstractVector) = T(map(myuparse, parameters))
_materialize_eos(T, parameters::AbstractDict) =
    T((myuparse(parameters[string(f)]) for f in fieldnames(T))...)

function materialize_press(config::Pressures)
    unit = myuparse(config.unit)
    return config.values .* unit
end

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

function materialize_dir(config::Directories, fixed::Union{Pressures,Volumes})
    return map(fixed.values) do value
        abspath(joinpath(expanduser(config.root), config.prefix * string(ustrip(value))))
    end
end
materialize_dir(config::EosFittingConfig) = materialize_dir(config.dirs, config.fixed)

function materialize end

end
