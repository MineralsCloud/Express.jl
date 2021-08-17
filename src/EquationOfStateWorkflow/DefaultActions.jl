module DefaultActions

using AbInitioSoftwareBase: save, load, extension
using AbInitioSoftwareBase.Inputs: Input, writetxt
using Dates: now, format
using EquationsOfStateOfSolids:
    EquationOfStateOfSolids, EnergyEquation, PressureEquation, Parameters, getparam
using EquationsOfStateOfSolids.Fitting: eosfit
using Logging: with_logger, current_logger
using Serialization: serialize, deserialize
using Unitful: ustrip, unit

using ...Express: Action, loadconfig
using ..EquationOfStateWorkflow: ScfOrOptim, Scf, CURRENT_CALCULATION

struct MakeCmd{T} <: Action{T} end

struct MakeInput{T} <: Action{T} end
function (x::MakeInput)(file, template::Input, args...)
    input = x(template, args...)
    mkpath(dirname(file))  # In case its parent directory is not created
    writetxt(file, input)
    return input
end

struct GetData{T} <: Action{T} end
function (::GetData{T})(outputs) where {T<:ScfOrOptim}
    raw = (parseoutput(T())(output) for output in outputs)  # `ntuple` cannot work with generators
    return collect(Iterators.filter(!isnothing, raw))  # A vector of pairs
end
function (x::GetData{T})(file, outputs) where {T}
    data = x(outputs)
    dict = Dict(
        "volume" => (ustrip ∘ first).(data),
        "energy" => (ustrip ∘ last).(data),
        "vunit" => string(unit(first(data).first)),
        "eunit" => string(unit(first(data).second)),
    )
    ext = lowercase(extension(file))
    if ext == "jls"
        open(file, "w") do io
            serialize(io, dict)
        end
    elseif ext in ("json", "yaml", "yml", "toml")
        save(file, dict)
    else
        error("unsupported file extension `$ext`!")
    end
end

function parseoutput end

struct FitEos{T} <: Action{T} end
function (x::FitEos{T})(
    data::AbstractVector{<:Pair},
    trial_eos::EnergyEquation,
) where {T<:ScfOrOptim}
    return eosfit(trial_eos, first.(data), last.(data))
end
function (x::FitEos{T})(outputs, trial_eos::EnergyEquation) where {T<:ScfOrOptim}
    data = GetData{T}()(outputs)
    if length(data) <= 5
        @info "pressures <= 5 may give unreliable results, run more if possible!"
    end
    return x(data, trial_eos)
end
function (x::FitEos{T})(cfgfile) where {T<:ScfOrOptim}
    config = loadconfig(cfgfile)
    outfiles = last.(iofiles(T(), cfgfile))
    saveto = joinpath(config.workdir, shortname(T) * "_eos.jls")
    trial_eos =
        T <: Scf ? config.trial_eos :
        deserialize(joinpath(config.workdir, shortname(Scf) * "_eos.jls"))
    eos = x(outfiles, EnergyEquation(trial_eos))
    SaveEos{T}()(saveto, eos)
    return eos
end

struct SaveEos{T} <: Action{T} end
function (::SaveEos{T})(file, eos::Parameters) where {T<:ScfOrOptim}
    ext = lowercase(extension(file))
    if ext == "jls"
        open(file, "w") do io
            serialize(io, eos)
        end
    elseif ext in ("json", "yaml", "yml", "toml")
        save(file, eos)
    else
        error("unsupported file extension `$ext`!")
    end
end
(x::SaveEos)(path, eos::EquationOfStateOfSolids) = x(path, getparam(eos))

struct LogMsg{T} <: Action{T} end
function (x::LogMsg{T})(start = true) where {T}
    startend = start ? "starts" : "ends"
    with_logger(current_logger()) do
        println("The calculation $T $startend at $(format(now(), "HH:MM:SS u dd, yyyy")).")
    end
end

end
