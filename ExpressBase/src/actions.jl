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

think(x::DownloadPotentials, template::Input) = Thunk(x, template)
think(x::DownloadPotentials, config::NamedTuple) = Thunk(x, config.template)

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

think(x::WriteInput, path) = Thunk(x, path)

struct RunCmd{T} <: Action{T}
    calculation::T
end

think(f::RunCmd, config::NamedTuple) = think(f, config.cli.mpi.np, config.files)
function think(x::RunCmd, np::Integer, files, kwargs...)
    jobsize = length(files)
    np = procs_per_job(np, jobsize)
    return map(files) do (input, output)
        Thunk(x, input, output; np=np, kwargs...)
    end
end
