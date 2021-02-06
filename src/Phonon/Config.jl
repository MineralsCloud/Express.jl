module Config

using Configurations: @option

using ...Express: myuparse

@option "template" struct DfptTemplate
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

@option "fit" struct PhononConfig
    templates::AbstractVector{DfptTemplate}
    fixed::Union{Pressures,Volumes}
    function PhononConfig(templates, fixed)
        @assert length(templates) >= 1
        if length(templates) != length(fixed.values)
            throw(
                DimensionMismatch(
                    "templates and pressures or volumes have different lengths!",
                ),
            )
        end
        return new(templates, fixed)
    end
end

function materialize_press_vol(config::Union{Pressures,Volumes})
    unit = myuparse(config.unit)
    return config.values .* unit
end

function checkconfig(config) end

function materialize end

end
