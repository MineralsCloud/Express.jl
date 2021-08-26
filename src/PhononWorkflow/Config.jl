module Config

using AbInitioSoftwareBase.Commands: CommandConfig
using Configurations: @option
using Unitful: ustrip

using ...Express: myuparse
using ...Config: Directories, @unit_vec_opt

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
    eos::String = "eos.jls"
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
    templates::Template
    fixed::Union{Pressures,Volumes}
    files::IOFiles = IOFiles()
    save::Save = Save()
    cli::CommandConfig
    function RuntimeConfig(templates, fixed, files, save, cli)
        if !isfile(template)
            @warn "I cannot find template file `$template`!"
        end
        return new(templates, fixed, files, save, cli)
    end
end

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
