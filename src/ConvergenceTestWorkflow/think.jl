using Thinkers: Thunk
using EasyConfig: Config as Conf

using .Config: StaticConfig, expand

import ExpressBase: think

think(action::CreateInput, conf::Conf) =
    collect(Thunk(action, conf.template, datum) for datum in conf.with)
think(action::ExtractData, conf::Conf) =
    collect(Thunk(action, file) for file in last.(conf.io))
think(action::SaveData, conf::Conf) = Thunk(action(conf.data.raw))
think(action::TestConvergence, conf::Conf) = Thunk(action, conf.threshold)
function think(action::Action{T}, config::StaticConfig) where {T}
    config = expand(config, T())
    return think(action, config::Conf)
end
