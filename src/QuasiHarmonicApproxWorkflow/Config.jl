module Config

using AbInitioSoftwareBase: save, load, parentdir
using Configurations: from_dict, @option
using Unitful: ustrip, @u_str
using ...Express: Action, myuparse

@option "pressures" struct Pressures
    values::AbstractVector
    unit::String
    function Pressures(values, unit = "GPa")
        if length(values) <= 5
            @info "less than 6 pressures may not fit accurately, consider adding more!"
        end
        return new(values, unit)
    end
end
Pressures(values::AbstractString, unit = "GPa") = Pressures(eval(Meta.parse(values)), unit)

@option "temperatures" struct Temperatures
    values::AbstractVector
    unit::String
    function Temperatures(values, unit = "K")
        @assert minimum(values) * myuparse(unit) >= 0u"K" "the minimum temperature is less than 0K!"
        return new(values, unit)
    end
end
Temperatures(values::AbstractString, unit = "K") =
    Temperatures(eval(Meta.parse(values)), unit)

@option struct Thermo
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

@option struct Directories
    root::String = pwd()
    pattern::String = "p=%.1f"
    group_by_step::Bool = false
    Directories(root, pattern, group_by_step) =
        new(abspath(expanduser(root)), pattern, group_by_step)
end

@option struct RuntimeConfig
    input::String
    inp_file_list::String
    static::String
    q_points::String
    temperatures::Temperatures
    pressures::Pressures
    thermo::Thermo = Thermo()
    dirs::Directories = Directories()
    calculation::String = "single"
    static_only::Bool = false
    order::UInt = 3
    energy_unit::String = "ry"
    function RuntimeConfig(
        input,
        inp_file_list,
        static,
        q_points,
        temperatures,
        pressures,
        thermo,
        dirs,
        calculation,
        static_only,
        order,
        energy_unit,
    )
        @assert lowercase(calculation) in
                ("single", "same phonon dos", "different phonon dos")
        @assert order in 3:5
        @assert lowercase(energy_unit) in ("ry", "ev")
        return new(
            input,
            inp_file_list,
            static,
            q_points,
            temperatures,
            pressures,
            thermo,
            dirs,
            lowercase(calculation),
            static_only,
            order,
            lowercase(energy_unit),
        )
    end
end

struct ExpandConfig{T} <: Action{T} end
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
    dict = Dict{String,Any}(
        "calculation" => config.calculation,
        "T_MIN" => ustrip(u"K", minimum(temperatures)),
        "NT" => length(temperatures),
        "DT" => ustrip(u"K", minimum(diff(temperatures))),
        "DT_SAMPLE" => ustrip(u"K", minimum(diff(temperatures))),
        "P_MIN" => ustrip(u"GPa", minimum(pressures)),
        "NTV" => length(pressures),
        "DELTA_P" => ustrip(u"GPa", minimum(diff(pressures))),
        "DELTA_P_SAMPLE" => ustrip(u"GPa", minimum(diff(pressures))),
        "input" => config.input,
        "thermodynamic_properties" => map(fieldnames(Thermo)) do f
            if getfield(config.thermo, f)
                string(f)
            end
        end,
        "energy_unit" => config.energy_unit,
        "high_verbosity" => true,
        "output_directory" => joinpath(config.dirs.root, "results"),
    )
    path = expanduser(joinpath(parentdir(dict["input"]), "settings.yaml"))
    save(path, dict)
    return (
        input = abspath(expanduser(dict["input"])),
        config = path,
        inp_file_list = config.inp_file_list,
        static = config.static,
        q_points = config.q_points,
    )
end

end
