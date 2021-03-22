module Config

using AbInitioSoftwareBase: save, load
using Configurations: @option

using ...Config: @vecunit

@vecunit Pressures "GPa" "pressures" begin
    if length(values) <= 5
        @info "less than 6 pressures may not fit accurately, consider adding more!"
    end
    if minimum(values) >= zero(eltype(values))
        @warn "for better fitting result, provide at least 1 negative pressure!"
    end
end

@vecunit Volumes "bohr^3" "volumes" begin
    if length(values) <= 5
        @info "less than 6 volumes may not fit accurately, consider adding more!"
    end
end

@vecunit Temperatures "K" "temperatures"

@vecunit SampledTemperatures "K"

@vecunit SampledPressures "GPa"

@option struct Sampled
    temperatures::SampledTemperatures
    pressures::SampledPressures
end

@option struct Directories
    root::String = pwd()
    prefix::String = "p="
    group_by_step::Bool = false
end

@option struct QhaConfig
    input::String
    temperatures::Temperatures
    pressures::Pressures
    sampled::Sampled
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
