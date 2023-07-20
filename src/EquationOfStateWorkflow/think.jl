using Thinkers: Thunk
using Unitful: dimension, @u_str

using .Config: StaticConfig, expand

import ExpressBase: think

think(action::CreateInput, config::NamedTuple) =
    collect(Thunk(action, config.template) for _ in Base.OneTo(length(config.fixed)))
think(action::WriteInput, config::NamedTuple) = think.(action, first.(config.io))
think(action::ExtractData, config::NamedTuple) =
    collect(Thunk(action, file) for file in last.(config.io))
think(action::SaveData, config::NamedTuple) = Thunk(action, config.data.raw)
think(action::SaveParameters, config::NamedTuple) = Thunk(action, config.save.parameters)
function think(action::FitEquationOfState, config::NamedTuple)
    trial_eos = if action.calculation isa SCF
        config.trial_eos
    else
        loadparameters(config.save.parameters)
    end
    return Thunk(action(EnergyEquation(trial_eos)), Set())
end
function think(action::Action, config::StaticConfig)
    config = expand(config, action.calculation)
    return think(action, config)
end
