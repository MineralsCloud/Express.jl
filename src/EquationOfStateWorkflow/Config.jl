module Config

using AbInitioSoftwareBase.Commands: CommandConfig
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
using Unitful: AbstractQuantity, ustrip

using ...Express: myuparse
using ...Config: Directories, @unit_vec_opt

@unit_vec_opt Pressures "GPa" "pressures" begin
    function (values, _)
        if length(values) <= 5
            @info "less than 6 pressures may not fit accurately, consider adding more!"
        end
        if minimum(values) >= zero(eltype(values))
            @warn "for better fitting result, provide at least 1 negative pressure!"
        end
    end
end

@unit_vec_opt Volumes "bohr^3" "volumes" begin
    function (values, _)
        if length(values) <= 5
            @info "less than 6 volumes may not fit accurately, consider adding more!"
        end
    end
end

@option struct TrialEquationOfState
    name::String
    parameters::Union{AbstractVector,AbstractDict}
end

@option struct NamingConvention
    input::String = "%s.in"
    output::String = "%s.out"
end

@option struct GeneratedFiles
    dirs::Directories = Directories()
    naming_convention::NamingConvention = NamingConvention()
end

@option struct RuntimeConfig
    template::String
    trial_eos::TrialEquationOfState
    fixed::Union{Pressures,Volumes,Nothing} = nothing
    files::GeneratedFiles = GeneratedFiles()
    recover::String = ""
    cli::CommandConfig
    function RuntimeConfig(template, trial_eos, fixed, files, recover, cli)
        if !isempty(recover)
            recover = abspath(expanduser(recover))
        end
        return new(template, trial_eos, fixed, files, recover, cli)
    end
end

function materialize(config::TrialEquationOfState)
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
