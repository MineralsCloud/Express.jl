using ExpressBase.Files: load
using AbInitioSoftwareBase: Input, getpseudodir, getpotentials
using Dates: format, now
using ExpressBase: Action, calculation
using Logging: with_logger, current_logger
using Pseudopotentials: download_potential
using Thinkers: Thunk

using ..Express: distribute_procs

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

think(x::DownloadPotentials, template::Input) = Thunk(x, template)
think(x::DownloadPotentials, config::NamedTuple) = Thunk(x, config.template)

struct RunCmd{T} <: Action{T} end

function think(f::RunCmd{T}, config::NamedTuple) where {T}
    return think(f, config.cli.mpi.np, config.files)
end
function think(x::RunCmd, np::Integer, files, kwargs...)
    jobsize = length(files)
    np = distribute_procs(np, jobsize)
    return map(files) do (input, output)
        Thunk(x, input, output; np=np, kwargs...)
    end
end
