using AbInitioSoftwareBase: save, load, extension
using AbInitioSoftwareBase.Inputs: Input, writetxt
using Dates: now, format
using EquationsOfStateOfSolids:
    EquationOfStateOfSolids, EnergyEquation, PressureEquation, Parameters, getparam
using EquationsOfStateOfSolids.Fitting: eosfit
using Logging: with_logger, current_logger
using Pseudopotentials: download_potential
using Serialization: serialize, deserialize
using SimpleWorkflows: AtomicJob
using Unitful: ustrip, unit

using ...Express: Action, calculation
using ..EquationOfStateWorkflow: ScfOrOptim, Scf, Optimization
using ..Shell: distprocs
using .Config: Volumes, ExpandConfig

struct DownloadPotentials{T} <: Action{T} end
function (x::DownloadPotentials)(template::Input)
    dir = getpseudodir(template)
    potentials = getpotentials(template)
    return map(potentials) do potential
        path = joinpath(dir, potential)
        if !isfile(path)
            download_potential(potential, path)
        end
    end
end

function buildjob(x::DownloadPotentials{T}, cfgfile) where {T}
    dict = load(cfgfile)
    config = ExpandConfig{T}()(dict)
    return AtomicJob(() -> x(config.template))
end

function getpseudodir end

function getpotentials end

struct MakeInput{T} <: Action{T} end
function (x::MakeInput)(file, template::Input, args...)
    input = x(template, args...)
    mkpath(dirname(file))  # In case its parent directory is not created
    writetxt(file, input)
    return input
end

function buildjob(x::MakeInput{Scf}, cfgfile)
    dict = load(cfgfile)
    config = ExpandConfig{Scf}()(dict)
    inputs = first.(config.files)
    if config.fixed isa Volumes
        return map(inputs, config.fixed) do input, volume
            AtomicJob(() -> x(input, config.template, volume, "Y-m-d_H:M:S"))
        end
    else  # Pressure
        return map(inputs, config.fixed) do input, pressure
            AtomicJob(
                function ()
                    trial_eos = PressureEquation(config.trial_eos)
                    return x(input, config.template, trial_eos, pressure, "Y-m-d_H:M:S")
                end,
            )
        end
    end
end
function buildjob(x::MakeInput{T}, cfgfile) where {T<:Optimization}
    dict = load(cfgfile)
    config = ExpandConfig{T}()(dict)
    inputs = first.(config.files)
    if config.fixed isa Volumes
        return map(inputs, config.fixed) do input, volume
            AtomicJob(() -> x(input, config.template, volume, "Y-m-d_H:M:S"))
        end
    else  # Pressure
        return map(inputs, config.fixed) do input, pressure
            AtomicJob(
                function ()
                    trial_eos = PressureEquation(
                        FitEos{Scf}()(
                            last.(ExpandConfig{Scf}()(dict).files),
                            EnergyEquation(config.trial_eos),
                        ),
                    )
                    return x(input, config.template, trial_eos, pressure, "Y-m-d_H:M:S")
                end,
            )
        end
    end
end

struct RunCmd{T} <: Action{T} end

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
    return collect(Iterators.filter(x -> x !== nothing, raw))  # A vector of pairs
end

function buildjob(x::GetData{T}, cfgfile) where {T}
    dict = load(cfgfile)
    config = ExpandConfig{T}()(dict)
    AtomicJob(
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

function buildjob(x::FitEos{T}, cfgfile) where {T<:ScfOrOptim}
    return AtomicJob(
        function ()
            dict = load(cfgfile)
            config = ExpandConfig{T}()(dict)
            outputs = last.(config.files)
            trial_eos =
                calculation(x) isa Scf ? config.trial_eos :
                deserialize(config.save_eos)[string(nameof(Scf))]
            eos = x(outputs, EnergyEquation(trial_eos))
            SaveEos{T}()(config.save_eos, eos)
            return eos
        end,
    )
end

struct SaveEos{T} <: Action{T} end
function (::SaveEos{T})(file, eos::Parameters) where {T}
    @assert extension(file) == "jls"
    dict = isfile(file) ? deserialize(file) : Dict()
    dict[string(nameof(T))] = eos
    serialize(file, dict)
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
