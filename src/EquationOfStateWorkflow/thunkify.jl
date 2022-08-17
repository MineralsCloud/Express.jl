thunkify(
    f::UpdateTemplate,
    template::Input,
    pressures::Pressures,
    trial_eos::PressureEquation,
) = map(pressure -> Thunk(f, template, trial_eos, pressure), pressures)
thunkify(f::UpdateTemplate, template::Input, volumes::Volumes) =
    map(volume -> Thunk(f, template, volume), volumes)
thunkify(f::UpdateTemplate{Scf}, config::NamedTuple) =
    thunkify(f, config.template, config.fixed, config.trial_eos)
function thunkify(f::UpdateTemplate{<:Optimization}, config::NamedTuple)
    if config.fixed isa Volumes
        thunkify(f, config.template, config.fixed)
    else  # Pressures
        trial_eos = PressureEquation(
            FitEos{Scf}()(last.(config.files), EnergyEquation(config.trial_eos)),
        )
        return thunkify(f, config.template, config.fixed, trial_eos)
    end
end
function thunkify(f::UpdateTemplate{T}, file::ConfigFile) where {T}
    raw_config = load(file)
    config = ExpandConfig{T}()(raw_config)
    return thunkify(f, config)
end
