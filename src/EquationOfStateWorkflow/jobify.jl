# MakeInput
function jobify(x::MakeInput, inputs, template::Input, pressures::Pressures, trial_eos)
    return map(inputs, pressures) do input, pressure
        Job(Thunk(x, input, template, trial_eos, pressure, "Y-m-d_H:M:S"))
    end
end
function jobify(x::MakeInput, inputs, template::Input, volumes::Volumes, args...)
    return map(inputs, volumes) do input, volume
        Job(Thunk(x, input, template, volume, "Y-m-d_H:M:S"))
    end
end
jobify(x::MakeInput{Scf}, config::NamedTuple) =
    jobify(x, first.(config.files), config.template, config.fixed, config.trial_eos)
function jobify(x::MakeInput{<:Optimization}, config::NamedTuple)
    trial_eos = PressureEquation(
        FitEos{Scf}()(last.(config.files), EnergyEquation(config.trial_eos)),
    )
    return jobify(x, first.(config.files), config.template, config.fixed, trial_eos)
end
function jobify(x::MakeInput{T}, file::ConfigFile) where {T}
    raw_config = load(file)
    config = ExpandConfig{T}()(raw_config)
    return jobify(x, config)
end
