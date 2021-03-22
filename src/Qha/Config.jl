    end
end

function checkconfig(config)
    for key in ("inp_file_list", "static", "q_points")
        @assert haskey(config, key) "`\"$key\"` was not found in config!"
    end
    list = load(config["inp_file_list"])["frequency_files"]
    _alert(list)
    for path in vcat(list, config["static"], config["q_points"])
        if !isfile(path)
            throw(SystemError("opening file \"$path\"", 2))
        end
    end
    @assert config["t_min"] >= 0 && config["nt"] >= 0 && config["npress"] >= 0
    if haskey(config, "use_shell") && haskey(config, "shell_args") && config["use_shell"]
        @assert config["shell_args"] isa AbstractDict
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
