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

function thunkify end

jobify(f::Action, args...) = Job.(thunkify(f, args...))

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

thunkify(x::DownloadPotentials, template::Input) = Thunk(x, template)
thunkify(x::DownloadPotentials, config::NamedTuple) = Thunk(x, config.template)
function thunkify(x::DownloadPotentials{T}, file::ConfigFile) where {T}
    raw_config = load(file)
    config = ExpandConfig{T}()(raw_config)
    return thunkify(x, config)
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

thunkify(x::LogTime) = Thunk(x, ())

struct RunCmd{T} <: Action{T} end

function thunkify(x::RunCmd, np::Integer, files, kwargs...)
    jobsize = length(files)
    np = distribute_procs(np, jobsize)
    return map(files) do (input, output)
        Thunk(x, input, output; np=np, kwargs...)
    end
end
