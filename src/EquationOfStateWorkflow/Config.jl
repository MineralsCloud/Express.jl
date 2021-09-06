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
using Formatting: sprintf1

using ...Express: Calculation, Action, myuparse

@option "pressures" struct Pressures
    values::AbstractVector
    unit::String
    function Pressures(values, unit = "GPa")
        if length(values) <= 5
            @info "less than 6 pressures may not fit accurately, consider adding more!"
        end
        if minimum(values) >= zero(eltype(values))
            @warn "for better fitting result, provide at least 1 negative pressure!"
        end
        return new(values, unit)
    end
end
Pressures(values::AbstractString, unit = "GPa") = Pressures(eval(Meta.parse(values)), unit)

@option "volumes" struct Volumes
    values::AbstractVector
    unit::String
    function Volumes(values, unit = "bohr^3")
        if length(values) <= 5
            @info "less than 6 volumes may not fit accurately, consider adding more!"
        end
        return new(values, unit)
    end
end
Volumes(values::AbstractString, unit = "bohr^3") = Volumes(eval(Meta.parse(values)), unit)

@option struct TrialEquationOfState
    type::String
    values::Union{AbstractVector,AbstractDict}
end

@option struct Directories
    root::String = pwd()
    pattern::String = "p=%.1f"
    group_by_step::Bool = false
    Directories(root, pattern, group_by_step) =
        new(abspath(expanduser(root)), pattern, group_by_step)
end

@option struct NamingPattern
    input::String = "%s.in"
    output::String = "%s.out"
end

@option struct Save
    raw::String = "raw.json"
    eos::String = "eos.jls"
    status::String = ""
end

@option struct IOFiles
    dirs::Directories = Directories()
    pattern::NamingPattern = NamingPattern()
end

@option struct RuntimeConfig
    template::String
    trial_eos::TrialEquationOfState
    fixed::Union{Pressures,Volumes}
    files::IOFiles = IOFiles()
    save::Save = Save()
    cli::CommandConfig
    function RuntimeConfig(template, trial_eos, fixed, files, save, cli)
        if !isfile(template)
            @warn "I cannot find template file `$template`!"
        end
        return new(template, trial_eos, fixed, files, save, cli)
    end
end

struct ExpandConfig{T} <: Action{T} end
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
    if trial_eos.values isa AbstractVector
        return T(map(myuparse, trial_eos.values)...)
    elseif trial_eos.values isa AbstractDict
        return T((myuparse(trial_eos.values[string(f)]) for f in propertynames(T))...)
    else
        @assert false "this is a bug!"
    end
end
function (::ExpandConfig)(fixed::Union{Pressures,Volumes})
    unit = myuparse(fixed.unit)
    return fixed.values .* unit
end
function (::ExpandConfig{T})(files::IOFiles, fixed::Union{Pressures,Volumes}) where {T}
    dirs = map(fixed.values) do value
        abspath(joinpath(files.dirs.root, sprintf1(files.dirs.pattern, value)))
    end
    return map(dirs) do dir
        type = string(nameof(T))
        in, out = sprintf1(files.pattern.input, type), sprintf1(files.pattern.output, type)
        joinpath(dir, in) => joinpath(dir, out)
    end
end
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
