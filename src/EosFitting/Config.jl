module Config

using AbInitioSoftwareBase.Cli: CliConfig
using Compat: isnothing
using Configurations: @option
using Crystallography: cellvolume
using EquationsOfStateOfSolids:
    Murnaghan1st,
    Murnaghan2nd,
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
using ...Config: @vecunit

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

@vecunit Pressures "GPa" "pressures" begin
    if length(values) <= 5
        @info "less than 6 pressures may not fit accurately, consider adding more!"
    end
    if minimum(values) >= zero(eltype(values))
        @warn "for better fitting result, provide at least 1 negative pressure!"
    end
end

@vecunit Volumes "bohr^3" "volumes" begin
    if length(values) <= 5
        @info "less than 6 volumes may not fit accurately, consider adding more!"
    end
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
    inv_opt::NumericalInversionOptions = NumericalInversionOptions()
    recover::String = ""
    cli::T
    function EosFittingConfig{T}(
        templates,
        trial_eos,
        fixed,
        dirs,
        inv_opt,
        recover,
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
        if !isempty(recover)
            recover = abspath(expanduser(recover))
        end
        return new(templates, trial_eos, fixed, dirs, inv_opt, recover, cli)
    end
end

function materialize(config::TrialEos)
    name = filter(c -> isletter(c) || isdigit(c), lowercase(config.name))
    T = if name in ("m", "murnaghan")
        Murnaghan1st
    elseif name == "m2" || occursin("murnaghan2", name)
        Murnaghan2nd
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
    return _materialize(T, config.parameters)
end
_materialize(T, parameters::AbstractVector) = T(map(myuparse, parameters)...)
_materialize(T, parameters::AbstractDict) =
    T((myuparse(parameters[string(f)]) for f in propertynames(T))...)
function materialize(config::Union{Pressures,Volumes})
    unit = myuparse(config.unit)
    return config.values .* unit
end
function materialize(config::Directories, fixed::Union{Pressures,Volumes})
    return map(fixed.values) do value
        abspath(joinpath(expanduser(config.root), config.prefix * string(ustrip(value))))
    end
end

end
