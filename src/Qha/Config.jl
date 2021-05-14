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
    sampled::Sampled
    thermo::Thermo
    dirs::Directories = Directories()
    static_only::Bool = false
    order::UInt = 3
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

function materialize(config)
    dict = Dict{String,Any}(
        "calculation" => "single",
        "T_MIN" => config["t_min"],
        "NT" => config["nt"],
        "DT" => config["dt"],
        "DT_SAMPLE" => config["dt_sample"],
        "P_MIN" => config["p_min"],
        "NTV" => config["npress"],
        "DELTA_P" => config["delta_p"],
        "DELTA_P_SAMPLE" => config["delta_p_sample"],
        "input" => config["input"],
        "thermodynamic_properties" => config["thermodynamic_properties"],
        "energy_unit" => lowercase(config["energy_unit"]),
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
