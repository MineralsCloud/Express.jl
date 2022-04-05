module Config

using AbInitioSoftwareBase: save, load, parentdir
using Configurations: from_dict, @option
using Unitful: ustrip, @u_str
using ...Express: Action, UnitfulVector, myuparse

@option struct Pressures <: UnitfulVector
    values::AbstractVector
    unit::String
    function Pressures(values, unit = "GPa")
        if length(values) <= 5
            @info "less than 6 pressures may not fit accurately, consider adding more!"
        end
        return new(values, unit)
    end
end

@option struct Temperatures <: UnitfulVector
    values::AbstractVector
    unit::String
    function Temperatures(values, unit = "K")
        @assert minimum(values) * myuparse(unit) >= 0u"K" "the minimum temperature is less than 0K!"
        return new(values, unit)
    end
end

@option mutable struct Thermo
    F::Bool = true
    G::Bool = true
    U::Bool = true
    H::Bool = true
    V::Bool = true
    Cp::Bool = true
    Cv::Bool = true
    alpha::Bool = true
    Bt::Bool = true
    Btp::Bool = true
    Bs::Bool = true
    gamma::Bool = true
end

# From https://github.com/mauro3/Parameters.jl/blob/ecbf8df/src/Parameters.jl#L554-L561
function Base.show(io::IO, x::Thermo)
    if get(io, :compact, false) || get(io, :typeinfo, nothing) == typeof(x)
        Base.show_default(IOContext(io, :limit => true), x)
    else
        # just dumping seems to give ok output, in particular for big data-sets:
        dump(IOContext(io, :limit => true), x, maxdepth = 1)
    end
end

@option struct Directories
    root::String = pwd()
    pattern::String = "p=%.1f"
    group_by_step::Bool = false
    Directories(root, pattern, group_by_step) =
        new(abspath(expanduser(root)), pattern, group_by_step)
end

@option struct RuntimeConfig
    recipe::String = "single qha"
    input::String
    inp_file_list::String
    static::String
    q_points::String
    temperatures::Temperatures
    pressures::Pressures
    thermo::Thermo = Thermo()
    dirs::Directories = Directories()
    static_only::Bool = false
    order::UInt = 3
    energy_unit::String = "ry"
    function RuntimeConfig(
        recipe,
        input,
        inp_file_list,
        static,
        q_points,
        temperatures,
        pressures,
        thermo,
        dirs,
        static_only,
        order,
        energy_unit,
    )
        @assert recipe in ("single qha", "multi qha")
        @assert order in 3:5
        @assert lowercase(energy_unit) in ("ry", "ev")
        return new(
            recipe,
            input,
            inp_file_list,
            static,
            q_points,
            temperatures,
            pressures,
            thermo,
            dirs,
            static_only,
            order,
            lowercase(energy_unit),
        )
    end
end

struct ExpandConfig{T} end
function (::ExpandConfig)(pressures::Pressures)
    unit = myuparse(pressures.unit)
    expanded = pressures.values .* unit
    if minimum(expanded) >= zero(eltype(expanded))  # values may have eltype `Any`
        @warn "for better fitting result, provide at least 1 negative pressure!"
    end
    return expanded
end
function (::ExpandConfig)(temperatures::Temperatures)
    unit = myuparse(temperatures.unit)
    return temperatures.values .* unit
end
function (x::ExpandConfig)(config::AbstractDict)
    config = from_dict(RuntimeConfig, config)
    temperatures = x(config.temperatures)
    pressures = x(config.pressures)
    if config.recipe == "single qha"
        calculation = "single"
    elseif config.recipe == "multi qha"
        calculation = "different phonon dos"
    else
        @assert false "this should never happen!"
    end
    dict = Dict{String,Any}(
        "calculation" => calculation,
        "T_MIN" => ustrip(u"K", minimum(temperatures)),
        "NT" => length(temperatures),
        "DT" => ustrip(u"K", minimum(diff(temperatures))),
        "DT_SAMPLE" => ustrip(u"K", minimum(diff(temperatures))),
        "P_MIN" => ustrip(u"GPa", minimum(pressures)),
        "NTV" => length(pressures),
        "DELTA_P" => ustrip(u"GPa", minimum(diff(pressures))),
        "DELTA_P_SAMPLE" => ustrip(u"GPa", minimum(diff(pressures))),
        "input" => config.input,
        "thermodynamic_properties" => collect(map(fieldnames(Thermo)) do f
            if getfield(config.thermo, f)
                string(f)
            end
        end),
        "energy_unit" => config.energy_unit,
        "high_verbosity" => true,
        "output_directory" => joinpath(config.dirs.root, "results"),
    )
    path = expanduser(joinpath(parentdir(dict["input"]), "settings.yaml"))
    save(path, dict)
    return (
        input = abspath(expanduser(dict["input"])),
        config = path,
        inp_file_list = abspath(expanduser(config.inp_file_list)),
        static = abspath(expanduser(config.static)),
        q_points = abspath(expanduser(config.q_points)),
    )
end

end
