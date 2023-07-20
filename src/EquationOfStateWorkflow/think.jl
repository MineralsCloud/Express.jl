using Thinkers: Thunk
using Unitful: dimension, @u_str

using .Config: StaticConfig, ExpandConfig

import ExpressBase: think

think(obj::CreateInput, config::NamedTuple) =
    collect(Thunk(obj, config.template) for _ in Base.OneTo(length(config.fixed)))
think(obj::WriteInput, config::NamedTuple) = think.(obj, first.(config.io))
think(obj::ExtractData, files::AbstractVector) = collect(Thunk(obj, file) for file in files)
think(obj::ExtractData, config::NamedTuple) = think(obj, last.(config.io))
think(obj::SaveData, config::NamedTuple) = Thunk(obj(config.save.raw), Set())
think(obj::SaveParameters, config::NamedTuple) = Thunk(obj(config.save.parameters), Set())
function think(obj::FitEquationOfState, config::NamedTuple)
    trial_eos =
        obj.calculation isa SCF ? config.trial_eos : loadparameters(config.save.parameters)
    return Thunk(obj(EnergyEquation(trial_eos)), Set())
end
function think(f::Action{T}, config::StaticConfig{T}) where {T}
    config = ExpandConfig(T())(config)
    return think(f, config)
end
