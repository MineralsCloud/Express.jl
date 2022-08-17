using AbInitioSoftwareBase: save, load, extension
using AbInitioSoftwareBase.Inputs: Input, writetxt
using EquationsOfStateOfSolids:
    EquationOfStateOfSolids, EnergyEquation, PressureEquation, Parameters, getparam
using EquationsOfStateOfSolids.Fitting: eosfit
using ExpressBase: Action, Scf, Optimization, calculation
import JLD2
using SimpleWorkflows.Jobs: Job
using SimpleWorkflows.Thunks: Thunk
using Unitful: Pressure, Volume, ustrip, unit

using ..Config: ConfigFile
using .Config: ExpandConfig, Pressures, Volumes

import Express: jobify

struct MakeInput{T} <: Action{T} end
(x::MakeInput)(file, template::Input, pressure::Pressure, eos::PressureEquation, args...) =
    x(file, x(template, pressure, eos, args...))
(x::MakeInput)(file, template::Input, volume::Volume, args...) =
    x(file, x(template, volume, args...))
function (x::MakeInput)(file, input::Input)
    mkpath(dirname(file))  # In case its parent directory is not created
    writetxt(file, input)
    return input
end

struct GetData{T} <: Action{T} end
function (x::GetData)(outputs)
    raw = (parseoutput(calculation(x))(output) for output in outputs)  # `ntuple` cannot work with generators
    return collect(Iterators.filter(x -> x !== nothing, raw))  # A vector of pairs
end

function jobify(x::GetData{T}, cfgfile) where {T}
    dict = load(cfgfile)
    config = ExpandConfig{T}()(dict)
    Job(
        function ()
            data = x(last.(config.files))
            savedata = Dict(
                "volume" => (ustrip ∘ first).(data),
                "energy" => (ustrip ∘ last).(data),
                "vunit" => string(unit(first(data).first)),
                "eunit" => string(unit(first(data).second)),
            )
            dict = isfile(config.save_raw) ? load(config.save_raw) : Dict()
            dict[string(nameof(T))] = savedata
            save(config.save_raw, dict)
            return data
        end,
    )
end

function parseoutput end

struct FitEos{T} <: Action{T} end
(x::FitEos)(data::AbstractVector{<:Pair}, trial_eos::EnergyEquation) =
    eosfit(trial_eos, first.(data), last.(data))
function (x::FitEos)(outputs, trial_eos::EnergyEquation)
    data = GetData{typeof(calculation(x))}()(outputs)
    if length(data) <= 5
        @info "pressures <= 5 may give unreliable results, run more if possible!"
    end
    return x(data, trial_eos)
end

function jobify(x::FitEos{T}, cfgfile) where {T<:Union{Scf,Optimization}}
    return Job(
        function ()
            dict = load(cfgfile)
            config = ExpandConfig{T}()(dict)
            outputs = last.(config.files)
            trial_eos =
                calculation(x) isa Scf ? config.trial_eos :
                JLD2.load(config.save_eos)[string(nameof(Scf))]
            eos = x(outputs, EnergyEquation(trial_eos))
            SaveEos{T}()(config.save_eos, eos)
            return eos
        end,
    )
end

struct SaveEos{T} <: Action{T} end
function (::SaveEos{T})(file, eos::Parameters) where {T}
    @assert extension(file) == "jld2"
    dict = isfile(file) ? JLD2.load(file) : Dict{String,Any}()
    dict[string(nameof(T))] = eos
    JLD2.save(file, dict)
end
(x::SaveEos)(path, eos::EquationOfStateOfSolids) = x(path, getparam(eos))
