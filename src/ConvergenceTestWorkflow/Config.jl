module Config

using AbInitioSoftwareBase.Commands: CommandConfig
using Configurations: from_dict, @option
using Formatting: sprintf1
using Unitful: @u_str

using ...Express: Calculation, Action, UnitfulVector

@option "ecutwfc" struct CutoffEnergies <: UnitfulVector
    values::AbstractVector
    unit::String = "Ry"
end

@option "k_mesh" struct MonkhorstPackGrids
    meshes::AbstractVector{<:AbstractVector{<:Integer}}
    shifts::AbstractVector{<:AbstractVector{<:Integer}} = fill([0, 0, 0], length(meshes))
    function MonkhorstPackGrids(meshes, shifts)
        if length(meshes) != length(shifts)
            throw(DimensionMismatch("`meshes` and `shifts` should have the same length!"))
        end
        for (mesh, shift) in zip(meshes, shifts)
            @assert all(mesh .>= 1)
            @assert all(0 .<= shift .<= 1)
        end
        return new(meshes, shifts)
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
    parameters::Union{CutoffEnergies,MonkhorstPackGrids}
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
function (::ExpandConfig)(energies::CutoffEnergies)
    unit = @u_str(energies.unit)
    return energies.values .* unit
end
function (::ExpandConfig)(x::MonkhorstPackGrids)
    return map(x.meshes, x.shifts) do mesh, shift
        (mesh, shift)
    end
end
function (::ExpandConfig{T})(files::IOFiles, energies::CutoffEnergies) where {T}
    dirs = map(energies.values) do value
        abspath(joinpath(files.dirs.root, sprintf1(files.dirs.pattern, value)))
    end
    return map(dirs) do dir
        type = string(nameof(T))
        in, out = sprintf1(files.pattern.input, type), sprintf1(files.pattern.output, type)
        joinpath(dir, in) => joinpath(dir, out)
    end
end
function (::ExpandConfig{T})(files::IOFiles, grids::MonkhorstPackGrids) where {T}
    values = zip(grids.meshes, grids.shifts)
    dirs = map(values) do value
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
