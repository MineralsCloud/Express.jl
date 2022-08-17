module Config

using AbInitioSoftwareBase.Commands: CommandConfig
using Configurations: from_dict, @option
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

using ...Config: DirStructure, iofiles, @sp

@sp Pressures "GPa" "pressures" begin
    function (values, unit)
        if length(values) <= 5
            @info "less than 6 pressures may not fit accurately, consider adding more!"
        end
        if minimum(values * unit) >= 0 * unit  # values may have eltype `Any`
            @warn "for better fitting result, provide at least 1 negative pressure!"
        end
    end
end

@sp Volumes "bohr^3" "volumes" begin
    function (values, _)
        if length(values) <= 5
            @info "less than 6 volumes may not fit accurately, consider adding more!"
        end
    end
end

@option struct TrialEquationOfState
    type::String
    values::AbstractVector
end

@option struct Save
    raw::String = "raw.json"
    eos::String = "eos.jld2"
    status::String = ""
end

@option struct RuntimeConfig
    recipe::String
    template::String
    trial_eos::TrialEquationOfState
    fixed::Union{Pressures,Volumes}
    dirstructure::DirStructure
    save::Save = Save()
    cli::CommandConfig
    function RuntimeConfig(recipe, template, trial_eos, fixed, dirstructure, save, cli)
        @assert recipe in ("eos",)
        if !isfile(template)
            @warn "I cannot find template file `$template`!"
        end
        return new(recipe, template, trial_eos, fixed, dirstructure, save, cli)
    end
end

struct ExpandConfig{T} end
function (::ExpandConfig)(trial_eos::TrialEquationOfState)
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
    return T(map(myuparse ∘ string, trial_eos.values)...)
end
(::ExpandConfig)(data::Union{Pressures,Volumes}) = collect(datum for datum in data)
(::ExpandConfig{T})(ds::DirStructure, fixed::Union{Pressures,Volumes}) where {T} =
    iofiles(ds, fixed.values, string(nameof(T)))
function (::ExpandConfig)(save::Save)
    return map((:raw, :eos, :status)) do f
        v = getfield(save, f)
        isempty(v) ? abspath(mktemp(; cleanup = false)[1]) : abspath(expanduser(v))
    end
end
function (x::ExpandConfig)(config::AbstractDict)
    config = from_dict(RuntimeConfig, config)
    save_raw, save_eos, save_status = x(config.save)
    return (
        template = x(config.template),
        trial_eos = x(config.trial_eos),
        fixed = x(config.fixed),
        root = config.files.dirs.root,
        files = x(config.files, config.fixed),
        save_raw = save_raw,
        save_eos = save_eos,
        save_status = save_status,
        cli = config.cli,
    )
end

end
