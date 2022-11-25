using ExpressBase.Files: load
using AbInitioSoftwareBase.Inputs: Input, getpseudodir, getpotentials
using Dates: format, now
using EasyJobsBase: DependentJob
using EasyJobsBase.Thunks: Thunk
using ExpressBase: Action, calculation
using Logging: with_logger, current_logger
using Pseudopotentials: download_potential

using ..Express: distribute_procs
using .Config: ConfigFile

function thunkify end

function jobify(f::Action, args...)
    thunks = thunkify(f, args...)
    if thunks isa AbstractArray
        return map(DependentJob, thunks)
    else
        return DependentJob(thunks)
    end
end

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

function thunkify(f::RunCmd{T}, config::NamedTuple) where {T}
    return thunkify(f, config.cli.mpi.np, config.files)
end
function thunkify(x::RunCmd, np::Integer, files, kwargs...)
    jobsize = length(files)
    np = distribute_procs(np, jobsize)
    return map(files) do (input, output)
        Thunk(x, input, output; np=np, kwargs...)
    end
end
