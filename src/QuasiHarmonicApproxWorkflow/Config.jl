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
    function Pressures(values, unit = "bohr^3")
        if length(values) <= 5
            @info "less than 6 volumes may not fit accurately, consider adding more!"
        end
        return new(values, unit)
    end
end
Volumes(values::AbstractString, unit = "bohr^3") = Volumes(eval(Meta.parse(values)), unit)

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
    temperatures::Temperatures
    pressures::Pressures
    thermo::Thermo
    dirs::Directories = Directories()
    calculation::String = "single"
    static_only::Bool = false
    order::UInt = 3
    energy_unit::String = "ry"
    function RuntimeConfig(
        input,
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
function (::ExpandConfig)(config::AbstractDict)
    config = from_dict(RuntimeConfig, config)
    dict = Dict{String,Any}(
        "calculation" => config.calculation,
        "T_MIN" => ustrip(u"K", minimum(config.temperatures)),
        "NT" => length(config.temperatures),
        "DT" => ustrip(u"K", minimum(diff(config.temperatures))),
        "DT_SAMPLE" => ustrip(u"K", minimum(diff(config.temperatures))),
        "P_MIN" => ustrip(u"GPa", minimum(config.pressures)),
        "NTV" => length(config.pressures),
        "DELTA_P" => ustrip(u"GPa", minimum(diff(config.pressures))),
        "DELTA_P_SAMPLE" => ustrip(u"GPa", minimum(diff(config.pressures))),
        "input" => config.input,
        "thermodynamic_properties" => map(fieldnames(config.thermo)) do f
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
        input = abspath(expanduser(config["input"])),
        config = path,
        inp_file_list = config["inp_file_list"],
        static = config["static"],
        q_points = config["q_points"],
    )
end

end
