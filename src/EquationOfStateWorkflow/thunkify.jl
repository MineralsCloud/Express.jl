# MakeInput
function thunkify(
    f::MakeInput,
    files,
    template,
    pressures::Pressures,
    eos::PressureEquation,
)
    return map(files, pressures) do file, pressure
        Thunk(f, file, template, pressure, eos)
    end
end
function thunkify(f::MakeInput, files, template, volumes::Volumes)
    return map(files, volumes) do file, volume
        Thunk(f, file, template, volume)
    end
end
thunkify(f::MakeInput{Scf}, config::NamedTuple) =
    thunkify(f, config.files, config.template, config.fixed, config.trial_eos)
function thunkify(f::MakeInput{<:Optimization}, config::NamedTuple)
    if config.fixed isa Volumes
        thunkify(f, config.files, config.template, config.fixed)
    else  # Pressures
        trial_eos = PressureEquation(
            FitEos{Scf}()(last.(config.files), EnergyEquation(config.trial_eos)),
        )
        return thunkify(f, config.files, config.template, config.fixed, trial_eos)
    end
end
# Action
function thunkify(f::Action{T}, file::ConfigFile) where {T}
    raw_config = load(file)
    config = ExpandConfig{T}()(raw_config)
    return thunkify(f, config)
end
