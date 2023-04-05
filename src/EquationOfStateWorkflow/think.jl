using Configurations: from_dict
using Thinkers: Thunk
using Unitful: dimension, @u_str

using ..Express: DownloadPotentials
using .Config: RuntimeConfig

import ..Express: think

function think(
    makeinput::MakeInput{Scf}, files, template::Input, pressures, eos::PressureEquation
)
    return map(files, pressures) do (input, _), pressure
        Thunk(makeinput, input, template, pressure, eos)
    end
end
function think(
    makeinput::MakeInput{<:Optimization},
    files,
    template::Input,
    pressures,
    eos::PressureEquation,
)
    return map(files, pressures) do (input, _), pressure
        Thunk(
            function (file, template, pressure, eos)
                eos = PressureEquation(
                    FitEquationOfState{Scf}()(last.(files), EnergyEquation(eos))
                )
                return makeinput(file, template, pressure, eos)
            end,
            input,
            template,
            pressure,
            eos,
        )
    end
end
function think(makeinput::MakeInput, files, template::Input, volumes)
    return map(files, volumes) do (input, _), volume
        Thunk(makeinput, input, template, volume)
    end
end
function think(makeinput::MakeInput, config::NamedTuple)
    if dimension(first(config.fixed)) == dimension(u"m^3")  # Volumes
        think(makeinput, config.files, config.template, config.fixed)
    else  # Pressures
        return think(
            makeinput, config.files, config.template, config.fixed, config.trial_eos
        )
    end
end
function think(save::SaveVolumeEnergy, config::NamedTuple)
    return Thunk(function ()
        data = save(last.(config.files))
        return save(config.save.ev, data)
    end)
end
function think(fit::FitEquationOfState, config::NamedTuple)
    return Thunk(function ()
        outputs = last.(config.files)
        trial_eos = if calculation(fit) isa Scf
            config.trial_eos
        else
            load(config.save.eos)
        end
        eos = fit(outputs, EnergyEquation(trial_eos))
        SaveParameters{typeof(calculation(fit))}()(config.save.eos, eos)
        return eos
    end)
end
function think(f::Action{T}, raw_config::RuntimeConfig) where {T}
    config = ExpandConfig{T}()(raw_config)
    return think(f, config)
end
function think(f::Action{T}, file::ConfigFile) where {T}
    dict = load(file)
    config = from_dict(RuntimeConfig, dict)
    return think(f, config)
end
