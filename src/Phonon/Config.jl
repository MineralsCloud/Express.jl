module Config

using AbInitioSoftwareBase.Cli: CliConfig
using Configurations: @option
using Unitful: ustrip

using ...Express: myuparse

@option struct DfptTemplate
    scf::String
    dfpt::String
    q2r::String
    disp::String
    function DfptTemplate(scf, dfpt, q2r, disp)
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

@option "volumes" struct Volumes
    values::AbstractVector
    unit::String = "bohr^3"
end

@option struct Directories
    root::String = pwd()
    prefix::String = "p="
    group_by_step::Bool = false
end

@option struct PhononConfig{T<:CliConfig}
    templates::AbstractVector{DfptTemplate}
    fixed::Union{Pressures,Volumes}
    dirs::Directories = Directories()
    recover::String = ""
    cli::T
    function PhononConfig{T}(templates, fixed, dirs, recover, cli::T) where {T}
        @assert length(templates) >= 1
        if length(templates) != 1
            if length(templates) != length(fixed.values)
                throw(
                    DimensionMismatch(
                        "templates and pressures or volumes have different lengths!",
                    ),
                )
            end
        end
        return new(templates, fixed, dirs, recover, cli)
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
