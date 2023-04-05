using AbInitioSoftwareBase.Inputs: Input
using EquationsOfStateOfSolids:
    EquationOfStateOfSolids, EnergyEquation, PressureEquation, Parameters, getparam
using EquationsOfStateOfSolids.Fitting: eosfit
using ExpressBase: Action, Scf, Optimization, calculation
using ExpressBase.Files: save, load, extension
using JLD2: JLD2
using Unitful: Pressure, Volume, ustrip, unit

using ..Config: ConfigFile
using .Config: ExpandConfig, Pressures, Volumes

struct MakeInput{T} <: Action{T} end
(makeinput::MakeInput)(
    path, template::Input, pressure::Pressure, eos::PressureEquation, args...
) = makeinput(path, makeinput(template, pressure, eos, args...))
(makeinput::MakeInput)(path, template::Input, volume::Volume, args...) =
    makeinput(path, makeinput(template, volume, args...))
function (makeinput::MakeInput)(path, input::Input)
    if isfile(path)
        @warn "File $path already exists! It will be overwritten!"
    end
    mkpath(dirname(path))  # In case its parent directory is not created
    open(path, "w") do io
        print(io, input)
    end
    return input
end

function readdata(x::Calculation, outputs)
    data = map(outputs) do output
        parseoutput(calculation(x), output)
    end
    return filter(!isnothing, data)
end

function parseoutput end

struct FitEquationOfState{T} <: Action{T} end
(fit::FitEquationOfState)(data::AbstractVector{<:Pair}, trial_eos::EnergyEquation) =
    eosfit(trial_eos, first.(data), last.(data))
function (fit::FitEquationOfState)(outputs, trial_eos::EnergyEquation)
    data = readdata(calculation(fit), outputs)
    if length(data) <= 5
        @info "pressures <= 5 may give unreliable results, run more if possible!"
    end
    return fit(data, trial_eos)
end

struct SaveEos{T} <: Action{T} end
function (::SaveEos{T})(file, eos::Parameters) where {T}
    @assert extension(file) == "jld2"
    dict = isfile(file) ? JLD2.load(file) : Dict{String,Any}()
    dict[string(nameof(T))] = eos
    return JLD2.save(file, dict)
end
(x::SaveEos)(path, eos::EquationOfStateOfSolids) = x(path, getparam(eos))

include("think.jl")
