using AbInitioSoftwareBase: Input, getpseudodir, listpotentials
using Pseudopotentials: download_potential
using Thinkers: Thunk

using ..Files: load

struct DownloadPotentials{T} <: Action{T}
    calculation::T
end
function (::DownloadPotentials)(template::Input)
    dir = getpseudodir(template)
    if !isdir(dir)
        mkpath(dir)
    end
    return map(listpotentials(template)) do potential
        path = joinpath(dir, potential)
        if !isfile(path)
            download_potential(potential, path)
        end
    end
end

think(x::DownloadPotentials, config) = Thunk(x, config.template)

struct WriteInput{T} <: Action{T}
    calculation::T
end
function (obj::WriteInput)(path, input::Input)
    if isfile(path)
        @warn "file `$path` already exists! It will be overwritten!"
    end
    if !isdir(dirname(path))
        @warn "parent directory of `$path` does not exist! It will be created!"
        mkpath(dirname(path))
    end
    open(path, "w") do io
        print(io, input)  # `print` is overridden by loading `*Format` packages.
    end
    return nothing
end
(obj::WriteInput)(path) = Base.Fix1(obj, path)

think(action::WriteInput, file) = Thunk(action, file)

struct RunCmd{T} <: Action{T}
    calculation::T
end

function think(x::RunCmd, config)
    njobs = length(config.io)
    np = procs_per_job(config.cli.mpi.np, njobs)
    return map(config.io) do (input, output)
        Thunk(x, input, output; np=np)
    end
end
