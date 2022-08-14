module Config

using AbInitioSoftwareBase.Commands: CommandConfig
using Configurations: from_dict, @option
using ExpressBase: Calculation, Action
using ExpressWorkflowMaker.Config: @vopt
using ExpressWorkflowMaker.Templates.Config: DirStructure, iofiles

@vopt CutoffEnergies "Ry" "ecutwfc"

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

@option struct Save
    raw::String = "raw.json"
    status::String = ""
end

@option struct RuntimeConfig
    recipe::String
    template::String
    parameters::Union{CutoffEnergies,MonkhorstPackGrids}
    dirstructure::DirStructure = DirStructure()
    save::Save = Save()
    cli::CommandConfig
    function RuntimeConfig(recipe, template, parameters, dirstructure, save, cli)
        @assert recipe in ("ecut", "k_mesh")
        if !isfile(template)
            @warn "I cannot find template file `$template`!"
        end
        return new(recipe, template, parameters, dirstructure, save, cli)
    end
end

struct ExpandConfig{T} end
(::ExpandConfig)(energies::CutoffEnergies) = energies.values .* energies.unit
function (::ExpandConfig)(x::MonkhorstPackGrids)
    return map(x.meshes, x.shifts) do mesh, shift
        (mesh, shift)
    end
end
(::ExpandConfig{T})(ds::DirStructure, energies::CutoffEnergies) where {T} =
    iofiles(ds, energies.values, string(nameof(T)))
(::ExpandConfig{T})(ds::DirStructure, grids::MonkhorstPackGrids) where {T} =
    iofiles(ds, zip(grids.meshes, grids.shifts), string(nameof(T)))
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
