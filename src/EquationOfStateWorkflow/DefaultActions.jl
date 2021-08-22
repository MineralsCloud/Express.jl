module DefaultActions

using AbInitioSoftwareBase: save, load, extension
using AbInitioSoftwareBase.Inputs: Input, writetxt
using Compat: isnothing
using Dates: now, format
using EquationsOfStateOfSolids:
    EquationOfStateOfSolids, EnergyEquation, PressureEquation, Parameters, getparam
using EquationsOfStateOfSolids.Fitting: eosfit
using Logging: with_logger, current_logger
using Serialization: serialize, deserialize
using SimpleWorkflows: AtomicJob
using Unitful: ustrip, unit

using ...Express: Action, calculation
using ..EquationOfStateWorkflow: ScfOrOptim, Scf
using ..Config: Volumes, ExpandConfig
using ...Shell: distprocs

struct RunCmd{T} <: Action{T} end

struct MakeInput{T} <: Action{T} end
function (x::MakeInput)(file, template::Input, args...)
    input = x(template, args...)
    mkpath(dirname(file))  # In case its parent directory is not created
    writetxt(file, input)
    return input
end

function buildjob(x::MakeInput{T}, cfgfile) where {T}
    dict = load(cfgfile)
    config = ExpandConfig{T}()(dict)
    inputs = first.(config.files)
    trial_eos = PressureEquation(
        calculation(x) isa Scf ? config.trial_eos :
        FitEos{Scf}()(last.(config.files), config.trial_eos),
    )
    if config.fixed isa Volumes
        return map(inputs, config.fixed) do input, volume
            AtomicJob(() -> x(input, config.template, volume, "Y-m-d_H:M:S"))
        end
    else  # Pressure
        return map(inputs, config.fixed) do input, pressure
            AtomicJob(() -> x(input, config.template, trial_eos, pressure, "Y-m-d_H:M:S"))
        end
    end
end

function buildjob(x::RunCmd{T}, cfgfile) where {T}
    dict = load(cfgfile)
    config = ExpandConfig{T}()(dict)
    np = distprocs(config.cli.mpi.np, length(config.files))
    return map(config.files) do (input, output)
        AtomicJob(() -> x(input, output; np = np))
    end
end

struct GetData{T} <: Action{T} end
function (x::GetData)(outputs)
    raw = (parseoutput(calculation(x))(output) for output in outputs)  # `ntuple` cannot work with generators
    return collect(Iterators.filter(!isnothing, raw))  # A vector of pairs
end
function (x::GetData)(file, outputs)
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
(x::FitEos)(data::AbstractVector{<:Pair}, trial_eos::EnergyEquation) =
    eosfit(trial_eos, first.(data), last.(data))
function (x::FitEos{T})(outputs, trial_eos::EnergyEquation) where {T<:ScfOrOptim}
    data = GetData{T}()(outputs)
    if length(data) <= 5
        @info "pressures <= 5 may give unreliable results, run more if possible!"
    end
    return x(data, trial_eos)
end

function buildjob(x::FitEos{T}, cfgfile) where {T<:ScfOrOptim}
    dict = load(cfgfile)
    config = ExpandConfig{T}()(dict)
    outputs = last.(config.files)
    saveto = joinpath(config.files.dirs.root, string(T) * "_eos.jls")
    trial_eos =
        calculation(x) isa Scf ? config.trial_eos :
        deserialize(joinpath(config.files.dirs.root, string(Scf) * "_eos.jls"))
    eos = x(outputs, EnergyEquation(trial_eos))
    SaveEos{T}()(saveto, eos)
    return eos
end

struct SaveEos{T} <: Action{T} end
function (::SaveEos)(file, eos::Parameters)
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
function (x::LogMsg)(; start = true)
    act = start ? "starts" : "ends"
    with_logger(current_logger()) do
        println(
            "The calculation $(calculation(x)) $act at $(format(now(), "HH:MM:SS u dd, yyyy")).",
        )
    end
end

end
