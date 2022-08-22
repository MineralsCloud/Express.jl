using Configurations: from_dict
using SimpleWorkflows.Thunks: Thunk
using Unitful: dimension, @u_str

using ..Express: DownloadPotentials
using .Config: RuntimeConfig

import ..Express: thunkify

# MakeInput
function thunkify(f::MakeInput, files, template::Input, pressures, eos::PressureEquation)
    return map(files, pressures) do file, pressure
        Thunk(f, file, template, pressure, eos)
    end
end
function thunkify(f::MakeInput, files, template::Input, volumes)
    return map(files, volumes) do file, volume
        Thunk(f, file, template, volume)
    end
end
function thunkify(f::MakeInput{Scf}, config::NamedTuple)
    return thunkify(f, config.files, config.template, config.fixed, config.trial_eos)
end
function thunkify(f::MakeInput{<:Optimization}, config::NamedTuple)
    if dimension(first(config.fixed)) == dimension(u"m^3")  # Volumes
        thunkify(f, config.files, config.template, config.fixed)
    else  # Pressures
        trial_eos = PressureEquation(
            FitEos{Scf}()(last.(config.files), EnergyEquation(config.trial_eos))
        )
        return thunkify(f, config.files, config.template, config.fixed, trial_eos)
    end
end
function thunkify(x::GetRawData{T}, config::NamedTuple) where {T}
    return Thunk(function ()
        data = x(last.(config.files))
        dict = isfile(config.save.ev) ? load(config.save.ev) : Dict()
        dict[string(nameof(T))] = data
        save(config.save.ev, dict)
        return data
    end, ())
end
# FitEos
function thunkify(x::FitEos, config::NamedTuple)
    return Thunk(function ()
        outputs = last.(config.files)
        trial_eos = if calculation(x) isa Scf
            config.trial_eos
        else
            JLD2.load(config.save.eos)[string(nameof(Scf))]
        end
        eos = x(outputs, EnergyEquation(trial_eos))
        SaveEos{typeof(calculation(x))}()(config.save.eos, eos)
        return eos
    end, ())
end
# Any Action
function thunkify(f::Action{T}, raw_config::RuntimeConfig) where {T}
    config = ExpandConfig{T}()(raw_config)
    return thunkify(f, config)
end
function thunkify(f::Action{T}, file::ConfigFile) where {T}
    dict = load(file)
    config = from_dict(RuntimeConfig, dict)
    return thunkify(f, config)
end
