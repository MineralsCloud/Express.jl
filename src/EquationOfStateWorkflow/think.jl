using Configurations: from_dict
using Thinkers: Thunk
using Unitful: dimension, @u_str

using .Config: RuntimeConfig

import ExpressBase: think

think(obj::CreateInput, template::Input, pressures, eos::PressureEquation) =
    collect(Thunk(obj, template, pressure, eos) for pressure in pressures)
think(obj::CreateInput, template::Input, volumes) =
    collect(Thunk(obj, template, volume) for volume in volumes)
function think(obj::CreateInput{SCF}, config::NamedTuple)
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
think(obj::ExtractData, files::AbstractVector) = collect(Thunk(obj, file) for file in files)
think(obj::ExtractData, config::NamedTuple) = think(obj, last.(config.files))
think(obj::SaveData, config::NamedTuple) = Thunk(obj(config.save.raw), Set())
think(obj::SaveParameters, config::NamedTuple) = Thunk(obj(config.save.parameters), Set())
function think(obj::FitEquationOfState, config::NamedTuple)
    trial_eos =
        obj.calculation isa SCF ? config.trial_eos : loadparameters(config.save.parameters)
    return Thunk(obj(EnergyEquation(trial_eos)), Set())
end
function think(f::Action{T}, raw_config::RuntimeConfig) where {T}
    config = ExpandConfig{T}()(raw_config)
    return think(f, config)
end
