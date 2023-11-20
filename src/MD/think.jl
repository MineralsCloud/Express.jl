using Thinkers: Thunk
using EasyConfig: Config as Conf
using Unitful: dimension, @u_str

using .Config: StaticConfig, expand

import ExpressBase: think

think(action::CreateInput, conf::Conf) =
    collect(Thunk(action(conf.template)) for _ in Base.OneTo(length(conf.at)))
think(action::CreateInput{<:MolecularDynamics}, conf::Conf) =
    collect(Thunk(action, template) for template in first.(conf.io))
think(action::ExtractCell, conf::Conf) =
    collect(Thunk(action, file) for file in last.(conf.io))
function think(action::Action{T}, config::StaticConfig) where {T}
    config = expand(config, T())
    return think(action, config::Conf)
end
