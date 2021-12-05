using AbInitioSoftwareBase: save, load, extension
using AbInitioSoftwareBase.Inputs: Input, writetxt
using Dates: now, format
using Logging: with_logger, current_logger
using Serialization: serialize, deserialize
using SimpleWorkflows: AtomicJob
using Unitful: ustrip, unit

using ...Express: Action, calculation
using ..Shell: distprocs
using .Config: ExpandConfig

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
    return map(inputs, config.parameters) do input, energy
        AtomicJob(() -> x(input, config.template, energy, "Y-m-d_H:M:S"))
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

function parseoutput end

struct GetData{T} <: Action{T} end
function (x::GetData)(outputs)
    raw = (parseoutput(output) for output in outputs)  # `ntuple` cannot work with generators
    return collect(Iterators.filter(x -> x !== nothing, raw))  # A vector of pairs
end

function buildjob(x::GetData{T}, cfgfile) where {T}
    dict = load(cfgfile)
    config = ExpandConfig{T}()(dict)
    return AtomicJob(
        function ()
            data = x(last.(config.files))
            saved = Dict("x" => (ustrip ∘ first).(data), "y" => (ustrip ∘ last).(data))
            save(config.save_raw, saved)
            return data
        end,
    )
end

struct TestConvergence{T} <: Action{T} end
(x::TestConvergence)(data) = isconvergent(data)

function buildjob(x::TestConvergence{T}, cfgfile) where {T}
    dict = load(cfgfile)
    config = ExpandConfig{T}()(dict)
    return AtomicJob(function ()
        data = GetData{T}()(last.(config.files))
        return x(data)
    end)
end

struct LogMsg{T} <: Action{T} end
function (x::LogMsg)(; start = true)
    act = start ? "starts" : "ends"
    with_logger(current_logger()) do
        println(
            "The calculation $(calculation(x)) $act at $(format(now(), "HH:MM:SS u dd, yyyy")).",
        )
    end
end
