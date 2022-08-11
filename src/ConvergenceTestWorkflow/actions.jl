using AbInitioSoftwareBase: save, load, extension
using AbInitioSoftwareBase.Inputs: Input, writetxt, getpseudodir, getpotentials
using Dates: now, format
using ExpressBase: Action, calculation
using Logging: with_logger, current_logger
using Pseudopotentials: download_potential
using SimpleWorkflows: Job
using Unitful: ustrip, unit

using ..Shell: distprocs
using .Config: ExpandConfig

struct DownloadPotentials{T} <: Action{T} end
function (x::DownloadPotentials)(template::Input)
    dir = getpseudodir(template)
    if !isdir(dir)
        mkpath(dir)
    end
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
    return Job(() -> x(config.template))
end

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
    if config.parameters isa AbstractVector{<:Tuple}
        return map(inputs, config.parameters) do input, (mesh, shift)
            Job(() -> x(input, config.template, mesh, shift, "Y-m-d_H:M:S"))
        end
    else
        return map(inputs, config.parameters) do input, energy
            Job(() -> x(input, config.template, energy, "Y-m-d_H:M:S"))
        end
    end
end

struct RunCmd{T} <: Action{T} end

function buildjob(x::RunCmd{T}, cfgfile) where {T}
    dict = load(cfgfile)
    config = ExpandConfig{T}()(dict)
    np = distprocs(config.cli.mpi.np, length(config.files))
    return map(config.files) do (input, output)
        Job(() -> x(input, output; np = np))
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
    return Job(function ()
        data = x(last.(config.files))
        saved = Dict("results" => (ustrip ∘ last).(data))
        save(config.save_raw, saved)
        return data
    end)
end

struct TestConvergence{T} <: Action{T} end
(x::TestConvergence)(data) = isconvergent(data)

function buildjob(x::TestConvergence{T}, cfgfile) where {T}
    dict = load(cfgfile)
    config = ExpandConfig{T}()(dict)
    return Job(function ()
        data = GetData{T}()(last.(config.files))
        saved = Dict("results" => (ustrip ∘ last).(data))
        save(config.save_raw, saved)
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
