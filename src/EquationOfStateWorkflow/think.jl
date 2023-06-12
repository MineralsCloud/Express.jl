using Configurations: from_dict
using Thinkers: Thunk
using Unitful: dimension, @u_str

using ..Express: DownloadPotentials
using .Config: RuntimeConfig

import ..Express: think

function think(
    makeinput::MakeInput, files, template::Input, pressures, eos::PressureEquation
)
    return map(files, pressures) do (input, _), pressure
        Thunk(makeinput, input, template, pressure, eos)
    end
end
function think(makeinput::MakeInput, files, template::Input, volumes)
    return map(files, volumes) do (input, _), volume
        Thunk(makeinput, input, template, volume)
    end
end
function think(makeinput::MakeInput{Scf}, config::NamedTuple)
    if dimension(first(config.fixed)) == dimension(u"m^3")  # Volumes
        return think(makeinput, config.files, config.template, config.fixed)
    else  # Pressures
        return think(
            makeinput, config.files, config.template, config.fixed, config.trial_eos
        )
    end
end
function think(makeinput::MakeInput{<:Optimization}, config::NamedTuple)
    if dimension(first(config.fixed)) == dimension(u"m^3")  # Volumes
        return think(makeinput, config.files, config.template, config.fixed)
    else  # Pressures
        return map(config.files, config.fixed) do (input, _), pressure
            Thunk(trial_eos -> makeinput(input, config.template, pressure, trial_eos))
        end
    end
end
think(save::SaveVolumeEnergy, config::NamedTuple) =
    Thunk(data -> save(config.save.ve, data))
function think(fit::FitEquationOfState, config::NamedTuple)
    return Thunk(function (data)
        trial_eos = if calculation(fit) isa Scf
            config.trial_eos
        else
            LoadParameters{Scf}()(config.save.eos)
        end
        return fit(data, EnergyEquation(trial_eos))
    end, Set())
end
function think(f::Action{T}, raw_config::RuntimeConfig) where {T}
    config = ExpandConfig{T}()(raw_config)
    return think(f, config)
end
