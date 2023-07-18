module Config

using Configurations: from_dict, @option
using ExpressBase: CommandConfig

using ...Config: SamplingPoints, DirStructure, iofiles

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

@option "pressures" struct Pressures <: SamplingPoints
    numbers::Vector{Float64}
    unit::FreeUnits
    Pressures(numbers, unit="GPa") = new(numbers, unit)
end

@option "volumes" struct Volumes <: SamplingPoints
    numbers::Vector{Float64}
    unit::FreeUnits
    Volumes(numbers, unit="bohr^3") = new(numbers, unit)
end

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
(::ExpandConfig)(data::Union{Pressures,Volumes}) = collect(datum for datum in data)
(::ExpandConfig{T})(ds::DirStructure, fixed::Union{Pressures,Volumes}) where {T} =
    iofiles(ds, fixed.numbers, string(nameof(T)))
function (::ExpandConfig)(save::Save)
    return map((:raw, :status)) do f
        v = getfield(save, f)
        isempty(v) ? abspath(mktemp(; cleanup=false)[1]) : abspath(expanduser(v))
    end
end
function (x::ExpandConfig)(config::AbstractDict)
    config = from_dict(RuntimeConfig, config)
    save_raw, save_status = x(config.save)
    return (
        template=x(config.template),
        fixed=x(config.fixed),
        root=config.files.dirs.root,
        files=x(config.files, config.fixed),
        save_raw=save_raw,
        save_status=save_status,
        cli=config.cli,
    )
end

end
