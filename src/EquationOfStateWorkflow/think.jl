using Configurations: from_dict
using Thinkers: Thunk
using Unitful: dimension, @u_str

using ..Express: DownloadPotentials
using .Config: RuntimeConfig

import ..Express: think

think(obj::CreateInput, template::Input, pressures, eos::PressureEquation) =
    collect(Thunk(obj, template, pressure, eos) for pressure in pressures)
think(obj::CreateInput, template::Input, volumes) =
    collect(Thunk(obj, template, volume) for volume in volumes)
function think(obj::CreateInput{Scf}, config::NamedTuple)
    if dimension(first(config.fixed)) == dimension(u"m^3")  # Volumes
        return think(obj, config.template, config.fixed)
    else  # Pressures
        return think(obj, config.template, config.fixed, config.trial_eos)
    end
end
function think(obj::CreateInput, config::NamedTuple)  # For optimizations
    if dimension(first(config.fixed)) == dimension(u"m^3")  # Volumes
        return think(obj, config.template, config.fixed)
    else  # Pressures
        return map(config.fixed) do pressure
            Thunk(trial_eos -> obj(config.template, pressure, trial_eos))
        end
    end
end
think(extract::ExtractData, files::AbstractVector) =
    map(input -> Thunk(extract, input), files)
think(extract::ExtractData, config::NamedTuple) = think(extract, last.(config.files))
think(save::SaveData, config::NamedTuple) = Thunk(data -> save(config.save.ev, data))
think(save::SaveParameters, config::NamedTuple) = Thunk(data -> save(config.save.eos, data))
function think(fit::FitEquationOfState, config::NamedTuple)
    return Thunk(function (data)
        trial_eos = if calculation(fit) isa Scf
            config.trial_eos
        else
            loadparameters(config.save.eos)
        end
        return fit(EnergyEquation(trial_eos))(data)
    end, Set())
end
function think(f::Action{T}, raw_config::RuntimeConfig) where {T}
    config = ExpandConfig{T}()(raw_config)
    return think(f, config)
end
