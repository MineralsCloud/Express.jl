module Config

using AbInitioSoftwareBase: save, load
using Configurations: from_dict, @option
using Unitful: ustrip, @u_str

using ...Config: Directories, @unit_vec_opt

@unit_vec_opt Pressures "GPa" "pressures" begin
    function (values, _)
        if length(values) <= 5
            @info "less than 6 pressures may not fit accurately, consider adding more!"
        end
        if minimum(values) >= zero(eltype(values))
            @warn "for better fitting result, provide at least 1 negative pressure!"
        end
    end
end

@unit_vec_opt Volumes "bohr^3" "volumes" begin
    function (values, _)
        if length(values) <= 5
            @info "less than 6 volumes may not fit accurately, consider adding more!"
        end
    end
end

@unit_vec_opt Temperatures "K" "temperatures" begin
    function (values, unit)
        @assert minimum(values) * unit >= 0u"K" "the minimum temperature is less than 0K!"
    end
end

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

@option struct QhaConfig
    input::String
    temperatures::Temperatures
    pressures::Pressures
    thermo::Thermo
    dirs::Directories = Directories()
    calculation::String = "single"
    static_only::Bool = false
    order::UInt = 3
    energy_unit::String = "ry"
    function QhaConfig(
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

function checkconfig(config)
    for key in ("inp_file_list", "static", "q_points")
        @assert haskey(config, key) "`\"$key\"` was not found in config!"
    end
    list = load(config["inp_file_list"])["frequency_files"]
    for path in vcat(list, config["static"], config["q_points"])
        if !isfile(path)
            throw(SystemError("opening file \"$path\"", 2))
        end
    end
    return
end

function materialize(config::AbstractDict)
    config = from_dict(QhaConfig, config)
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
        "output_directory" => joinpath(config["workdir"], "results"),
    )
    path = expanduser(joinpath(dirname(config["input"]), "settings.yaml"))
    save(path, dict)
    return (
        input = expanduser(config["input"]),
        config = path,
        inp_file_list = config["inp_file_list"],
        static = config["static"],
        q_points = config["q_points"],
    )
end

end
