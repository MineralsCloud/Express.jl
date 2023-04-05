using Configurations: from_dict
using Thinkers: Thunk
using Unitful: dimension, @u_str

using ..Express: DownloadPotentials
using .Config: RuntimeConfig

import ..Express: think

# MakeInput
function think(f::MakeInput{Scf}, files, template::Input, pressures, eos::PressureEquation)
    return map(files, pressures) do (input, _), pressure
        Thunk(f, input, template, pressure, eos)
    end
end
function think(
    f::MakeInput{<:Optimization}, files, template::Input, pressures, eos::PressureEquation
)
    return map(files, pressures) do (input, _), pressure
        Thunk(
            function (file, template, pressure, eos)
                eos = PressureEquation(
                    FitEquationOfState{Scf}()(last.(files), EnergyEquation(eos))
                )
                return f(file, template, pressure, eos)
            end,
            input,
            template,
            pressure,
            eos,
        )
    end
end
function think(f::MakeInput, files, template::Input, volumes)
    return map(files, volumes) do (input, _), volume
        Thunk(f, input, template, volume)
    end
end
function think(f::MakeInput, config::NamedTuple)
    if dimension(first(config.fixed)) == dimension(u"m^3")  # Volumes
        think(f, config.files, config.template, config.fixed)
    else  # Pressures
        return think(f, config.files, config.template, config.fixed, config.trial_eos)
    end
end
function think(x::GetRawData{Scf}, config::NamedTuple)
    return Thunk(function ()
        @show config.files
        data = x(last.(config.files))
        dict = isfile(config.save.ev) ? load(config.save.ev) : Dict()
        dict[string(nameof(Scf))] = data
        save(config.save.ev, dict)
        return data
    end)
end
function think(x::GetRawData{T}, config::NamedTuple) where {T<:Optimization}
    return Thunk(function ()
        data = x(last.(config.scffiles))
        dict = isfile(config.save.ev) ? load(config.save.ev) : Dict()
        dict[string(nameof(T))] = data
        save(config.save.ev, dict)
        return data
    end)
end
# FitEos
function think(x::FitEquationOfState, config::NamedTuple)
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
    end)
end
# Any Action
function think(f::Action{T}, raw_config::RuntimeConfig) where {T}
    config = ExpandConfig{T}()(raw_config)
    return think(f, config)
end
function think(f::Action{T}, file::ConfigFile) where {T}
    dict = load(file)
    config = from_dict(RuntimeConfig, dict)
    return think(f, config)
end
