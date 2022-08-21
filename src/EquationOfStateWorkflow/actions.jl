using AbInitioSoftwareBase: save, load, extension
using AbInitioSoftwareBase.Inputs: Input, writetxt
using EquationsOfStateOfSolids:
    EquationOfStateOfSolids, EnergyEquation, PressureEquation, Parameters, getparam
using EquationsOfStateOfSolids.Fitting: eosfit
using ExpressBase: Action, Scf, Optimization, calculation
using JLD2: JLD2
using Unitful: Pressure, Volume, ustrip, unit

using ..Config: ConfigFile
using .Config: ExpandConfig, Pressures, Volumes

struct MakeInput{T} <: Action{T} end
function (x::MakeInput)(
    file, template::Input, pressure::Pressure, eos::PressureEquation, args...
)
    return x(file, x(template, pressure, eos, args...))
end
function (x::MakeInput)(file, template::Input, volume::Volume, args...)
    return x(file, x(template, volume, args...))
end
function (x::MakeInput)(file, input::Input)
    mkpath(dirname(file))  # In case its parent directory is not created
    writetxt(file, input)
    return input
end

struct GetRawData{T} <: Action{T} end
function (x::GetRawData)(outputs)
    raw = map(outputs) do output
        parseoutput(calculation(x))(output)
    end
    return filter(!isnothing, raw)
end

function parseoutput end

struct FitEos{T} <: Action{T} end
function (x::FitEos)(data::AbstractVector{<:Pair}, trial_eos::EnergyEquation)
    return eosfit(trial_eos, first.(data), last.(data))
end
function (x::FitEos)(outputs, trial_eos::EnergyEquation)
    data = GetRawData{typeof(calculation(x))}()(outputs)
    if length(data) <= 5
        @info "pressures <= 5 may give unreliable results, run more if possible!"
    end
    return x(data, trial_eos)
end

struct SaveEos{T} <: Action{T} end
function (::SaveEos{T})(file, eos::Parameters) where {T}
    @assert extension(file) == "jld2"
    dict = isfile(file) ? JLD2.load(file) : Dict{String,Any}()
    dict[string(nameof(T))] = eos
    return JLD2.save(file, dict)
end
(x::SaveEos)(path, eos::EquationOfStateOfSolids) = x(path, getparam(eos))

include("thunkify.jl")
