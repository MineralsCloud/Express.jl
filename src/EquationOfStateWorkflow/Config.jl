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
using Formatting: sprintf1
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
    type::String
    values::Union{AbstractVector,AbstractDict}
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

function materialize(trial_eos::TrialEquationOfState)
    type = filter(c -> isletter(c) || isdigit(c), lowercase(trial_eos.type))
    T = if type in ("m", "murnaghan")
        Murnaghan1st
    elseif type == "m2" || occursin("murnaghan2", type)
        Murnaghan2nd
    elseif type == "bm2" || occursin("birchmurnaghan2", type)
        BirchMurnaghan2nd
    elseif type == "bm3" || occursin("birchmurnaghan3", type)
        BirchMurnaghan3rd
    elseif type == "bm4" || occursin("birchmurnaghan4", type)
        BirchMurnaghan4th
    elseif type == "pt2" || occursin("poiriertarantola2", type)
        PoirierTarantola2nd
    elseif type == "pt3" || occursin("poiriertarantola3", type)
        PoirierTarantola3rd
    elseif type == "pt4" || occursin("poiriertarantola4", type)
        PoirierTarantola4th
    elseif type == "v" || occursin("vinet", type)
        Vinet
    else
        error("unsupported eos name `\"$type\"`!")
    end
    if trial_eos.values isa AbstractVector
        return T(map(myuparse, trial_eos.values)...)
    elseif trial_eos.values isa AbstractDict
        return T((myuparse(trial_eos.values[string(f)]) for f in propertynames(T))...)
    else
        @assert false "this is a bug!"
    end
end
function materialize(fixed::Union{Pressures,Volumes})
    unit = myuparse(fixed.unit)
    return fixed.values .* unit
end
function materialize(dirs::Directories, fixed::Union{Pressures,Volumes})
    return map(fixed.values) do value
        abspath(joinpath(dirs.root, sprintf1(dirs.naming_convention, ustrip(value))))
    end
end

end
