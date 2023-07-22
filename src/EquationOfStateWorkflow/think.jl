using Thinkers: Thunk
using EasyConfig: Config as Conf
using Unitful: dimension, @u_str

using .Config: StaticConfig, expand

import ExpressBase: think

function think(action::ComputeVolume, conf::Conf)
    trial_eos = if action.calculation isa SelfConsistentField
        conf.trial_eos
    else
        loadparameters(conf.data.eos_params)
    end
    return collect(Thunk(action, pressure, trial_eos) for pressure in conf.at)
end
think(action::CreateInput, conf::Conf) =
    collect(Thunk(action, conf.template) for _ in Base.OneTo(length(conf.at)))
think(action::WriteInput, conf::Conf) =
    collect(think(action, file) for file in first.(conf.io))
think(action::ExtractData, conf::Conf) =
    collect(Thunk(action, file) for file in last.(conf.io))
think(action::SaveData, conf::Conf) = Thunk(action, conf.data.raw)
think(action::SaveParameters, conf::Conf) = Thunk(action, conf.data.eos_params)
function think(action::FitEquationOfState, conf::Conf)
    function fit(data)
        trial_eos = if action.calculation isa SelfConsistentField
            conf.trial_eos
        else
            loadparameters(conf.data.eos_params)
        end
        eos = EnergyEquation(trial_eos)
        return action(eos, data)
    end
    return Thunk(fit, Set())
end
function think(action::Action{T}, config::StaticConfig) where {T}
    config = expand(config, T())
    return think(action, config::Conf)
end
