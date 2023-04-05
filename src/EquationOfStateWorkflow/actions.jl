using AbInitioSoftwareBase.Inputs: Input
using EquationsOfStateOfSolids:
    EquationOfStateOfSolids, EnergyEquation, PressureEquation, Parameters, getparam
using EquationsOfStateOfSolids.Fitting: eosfit
using ExpressBase: Action, Scf, Optimization, calculation
using ExpressBase.Files: save, load, extension
using JLD2: JLD2
using Unitful: Pressure, Volume, ustrip, unit

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

struct SaveVolumeEnergy{T} <: Action{T} end
function (::SaveVolumeEnergy{T})(path, outputs) where {T}
    data = readdata(T(), outputs)
    dict = Dict("volume" => first.(data), "energy" => last.(data))
    return save(path, dict)
end

struct FitEquationOfState{T} <: Action{T} end
(fit::FitEquationOfState)(data::AbstractVector{<:Pair}, trial_eos::EnergyEquation) =
    eosfit(trial_eos, first.(data), last.(data))
function (fit::FitEquationOfState{T})(outputs, trial_eos::EnergyEquation) where {T}
    data = readdata(T(), outputs)
    if length(data) <= 5
        @info "pressures <= 5 may give unreliable results, run more if possible!"
    end
    return fit(data, trial_eos)
end

struct SaveParameters{T} <: Action{T} end
function (::SaveParameters)(path, parameters::Parameters)
    dict = Dict(
        "type" => string(typeof(parameters)),
        "params" => collect(
            string(getproperty(parameters, name)) for name in propertynames(parameters)
        ),
    )
    return save(path, dict)
end
(x::SaveParameters)(path, eos::EquationOfStateOfSolids) = x(path, getparam(eos))

include("think.jl")
