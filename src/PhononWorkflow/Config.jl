module Config

using AbInitioSoftwareBase.Commands: CommandConfig
using Configurations: from_dict, @option
using ExpressBase: Scf, LatticeDynamics, Dfpt, RealSpaceForceConstants, Action
using ExpressWorkflowMaker.Templates.Config: DirStructure, iofiles
using Formatting: sprintf1
using Unitful: ustrip

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

@vopt Pressures "GPa" "pressures"

@vopt Volumes "bohr^3" "volumes"

@option struct Save
    raw::String = "raw.json"
    status::String = ""
end

@option struct RuntimeConfig
    recipe::String
    template::Template
    fixed::Union{Pressures,Volumes}
    dirstructure::DirStructure = DirStructure()
    save::Save = Save()
    cli::CommandConfig
    function RuntimeConfig(recipe, template, fixed, dirstructure, save, cli)
        @assert recipe in ("phonon dispersion", "vdos")
        for i in 1:nfields(template)
            if !isfile(getfield(template, i))
                @warn "I cannot find template file `$template`!"
            end
        end
        return new(recipe, template, fixed, dirstructure, save, cli)
    end
end

struct ExpandConfig{T} end
(::ExpandConfig)(fixed::Union{Pressures,Volumes}) = fixed.values .* fixed.unit
function (::ExpandConfig)(save::Save)
    return map((:raw, :status)) do f
        v = getfield(save, f)
        isempty(v) ? abspath(mktemp(; cleanup = false)[1]) : abspath(expanduser(v))
    end
end
(::ExpandConfig{T})(ds::DirStructure, fixed::Union{Pressures,Volumes}) where {T} =
    iofiles(ds, fixed.values, string(nameof(T)))
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
