using AbInitioSoftwareBase: load
using AbInitioSoftwareBase.Inputs: Input, getpseudodir, getpotentials
using Dates: format, now
using ExpressBase: Action, calculation
using Logging: with_logger, current_logger
using Pseudopotentials: download_potential
using SimpleWorkflows.Jobs: Job
using SimpleWorkflows.Thunks: Thunk

using ..Express: distribute_procs
using .Config: ConfigFile

struct DownloadPotentials{T} <: Action{T} end
function (::DownloadPotentials)(template::Input)
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

jobify(x::DownloadPotentials, template::Input) = Job(Thunk(x, template))
jobify(x::DownloadPotentials, config::NamedTuple) = Job(Thunk(x, config.template))
function jobify(x::DownloadPotentials{T}, file::ConfigFile) where {T}
    raw_config = load(file)
    config = ExpandConfig{T}()(raw_config)
    return jobify(x, config)
end

struct LogTime{T} <: Action{T} end
function (x::LogTime)()
    with_logger(current_logger()) do
        println(
            "The calculation ",
            calculation(x),
            "starts at ",
            format(now(), "HH:MM:SS u dd, yyyy"),
            '.',
        )
    end
end

jobify(x::LogTime) = Job(Thunk(x, ()))

struct RunCmd{T} <: Action{T} end

function jobify(x::RunCmd, np::Integer, files, kwargs...)
    jobsize = length(files)
    np = distribute_procs(np, jobsize)
    return map(files) do (input, output)
        Job(Thunk(x, input, output; np = np, kwargs...))
    end
end
