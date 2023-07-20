using Thinkers: Thunk
using Unitful: dimension, @u_str

using .Config: StaticConfig, expand

import ExpressBase: think

think(obj::CreateInput, config::NamedTuple) =
    collect(Thunk(obj, config.template) for _ in Base.OneTo(length(config.fixed)))
think(obj::WriteInput, config::NamedTuple) = think.(obj, first.(config.io))
think(obj::ExtractData, config::NamedTuple) =
    collect(Thunk(obj, file) for file in last.(config.io))
think(obj::SaveData, config::NamedTuple) = Thunk(obj, config.data.raw)
think(obj::SaveParameters, config::NamedTuple) = Thunk(obj, config.save.parameters)
function think(obj::FitEquationOfState, config::NamedTuple)
    trial_eos =
        obj.calculation isa SCF ? config.trial_eos : loadparameters(config.save.parameters)
    return Thunk(obj(EnergyEquation(trial_eos)), Set())
end
function think(action::Action, config::StaticConfig)
    config = expand(config, action.calculation)
    return think(action, config)
end
