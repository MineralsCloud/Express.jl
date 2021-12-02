module Config

using AbInitioSoftwareBase.Commands: CommandConfig
using Configurations: from_dict, @option
using Formatting: sprintf1

using ...Express: Calculation, Action, UnitfulVector, myuparse

@option "ecutwfc" struct CutoffEnergy <: UnitfulVector
    values::AbstractVector
    unit::String = "Ry"
end

@option "k_mesh" struct KMesh <: UnitfulVector
    mesh::AbstractVector
    is_shift::AbstractVector = [0, 0, 0]
    function KMesh(mesh, is_shift)
        @assert all(mesh .>= 1)
        @assert all(0 <= x <= 1 for x in is_shift)
        return new(mesh, is_shift)
    end
end

@option struct Directories
    root::String = pwd()
    pattern::String = "e=%.1f"
end

@option struct FileNamePatterns
    input::String = "%s.in"
    output::String = "%s.out"
end

@option struct Save
    raw::String = "raw.json"
    status::String = ""
end

@option struct IOFiles
    dirs::Directories = Directories()
    pattern::FileNamePatterns = FileNamePatterns()
end

@option struct RuntimeConfig
    recipe::String
    template::String
    parameters::Union{CutoffEnergy,KMesh}
    files::IOFiles = IOFiles()
    save::Save = Save()
    cli::CommandConfig
    function RuntimeConfig(recipe, template, parameters, files, save, cli)
        @assert recipe in ("ecut", "k_mesh")
        if !isfile(template)
            @warn "I cannot find template file `$template`!"
        end
        return new(recipe, template, parameters, files, save, cli)
    end
end

struct ExpandConfig{T} <: Action{T} end
function (::ExpandConfig{T})(
    files::IOFiles,
    parameters::Union{CutoffEnergy,KMesh},
) where {T}
    dirs = map(parameters.values) do value
        abspath(joinpath(files.dirs.root, sprintf1(files.dirs.pattern, value)))
    end
    return map(dirs) do dir
        type = string(nameof(T))
        in, out = sprintf1(files.pattern.input, type), sprintf1(files.pattern.output, type)
        joinpath(dir, in) => joinpath(dir, out)
    end
end
function (::ExpandConfig)(save::Save)
    return map((:raw, :status)) do f
        v = getfield(save, f)
        isempty(v) ? abspath(mktemp(; cleanup = false)[1]) : abspath(expanduser(v))
    end
end
function (x::ExpandConfig)(config::AbstractDict)
    config = from_dict(RuntimeConfig, config)
    save_raw, save_status = x(config.save)
    return (
        template = x(config.template),
        parameters = x(config.parameters),
        root = config.files.dirs.root,
        files = x(config.files, config.parameters),
        save_raw = save_raw,
        save_status = save_status,
        cli = config.cli,
    )
end

end
