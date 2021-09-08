module Config

using AbInitioSoftwareBase.Commands: CommandConfig
using Configurations: from_dict, @option
using Formatting: sprintf1
using Unitful: ustrip
using ...Express: Action, myuparse
using ..PhononWorkflow: Scf, Dfpt, RealSpaceForceConstants, LatticeDynamics

@option struct Template
    scf::String
    dfpt::String
    q2r::String
    disp::String
    function Template(scf, dfpt, q2r, disp)
        if !isfile(scf)
            @warn "file \"$scf\" is not reachable, be careful!"
        end
        if !isfile(dfpt)
            @warn "file \"$dfpt\" is not reachable, be careful!"
        end
        if !isfile(q2r)
            @warn "file \"$q2r\" is not reachable, be careful!"
        end
        if !isfile(disp)
            @warn "file \"$disp\" is not reachable, be careful!"
        end
        return new(scf, dfpt, q2r, disp)
    end
end

@option "pressures" struct Pressures
    values::AbstractVector
    unit::String = "GPa"
end
Pressures(values::AbstractString, unit = "GPa") = Pressures(eval(Meta.parse(values)), unit)

@option "volumes" struct Volumes
    values::AbstractVector
    unit::String = "bohr^3"
end
Volumes(values::AbstractString, unit = "bohr^3") = Volumes(eval(Meta.parse(values)), unit)

@option struct Directories
    root::String = pwd()
    pattern::String = "p=%.1f"
    group_by_step::Bool = false
    Directories(root, pattern, group_by_step) =
        new(abspath(expanduser(root)), pattern, group_by_step)
end

@option struct Save
    raw::String = "raw.json"
    status::String = ""
end

@option struct NamingPattern
    input::String = "%s.in"
    output::String = "%s.out"
end

@option struct IOFiles
    dirs::Directories = Directories()
    pattern::NamingPattern = NamingPattern()
end

@option struct RuntimeConfig
    template::Template
    fixed::Union{Pressures,Volumes}
    files::IOFiles = IOFiles()
    save::Save = Save()
    cli::CommandConfig
    function RuntimeConfig(template, fixed, files, save, cli)
        for i in 1:nfields(template)
            if !isfile(getfield(template, i))
                @warn "I cannot find template file `$template`!"
            end
        end
        return new(template, fixed, files, save, cli)
    end
end

struct ExpandConfig{T} <: Action{T} end
function (::ExpandConfig)(fixed::Union{Pressures,Volumes})
    unit = myuparse(fixed.unit)
    return fixed.values .* unit
end
function (::ExpandConfig)(save::Save)
    return map((:raw, :status)) do f
        v = getfield(save, f)
        isempty(v) ? abspath(mktemp(; cleanup = false)[1]) : abspath(expanduser(v))
    end
end
function (x::ExpandConfig)(files::IOFiles, fixed::Union{Pressures,Volumes})
    dirs = map(fixed.values) do value
        abspath(joinpath(files.dirs.root, sprintf1(files.dirs.pattern, value)))
    end
    return map((Scf, Dfpt, RealSpaceForceConstants, LatticeDynamics)) do type
        map(dirs) do dir
            in, out = sprintf1(files.pattern.input, string(nameof(type))),
            sprintf1(files.pattern.output, string(nameof(type)))
            joinpath(dir, in) => joinpath(dir, out)
        end
    end
end
function (x::ExpandConfig)(config::AbstractDict)
    config = from_dict(RuntimeConfig, config)
    save_raw, save_status = x(config.save)
    return (
        template = x(config.template),
        fixed = x(config.fixed),
        root = config.files.dirs.root,
        files = x(config.files, config.fixed),
        save_raw = save_raw,
        save_status = save_status,
        cli = config.cli,
    )
end

end
