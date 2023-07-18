using AbInitioSoftwareBase: Input, getpseudodir, listpotentials
using Pseudopotentials: download_potential
using Thinkers: Thunk

using ..ExpressBase: distribute_procs
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
    path::String
end
function (obj::WriteInput)(input::Input)
    path = obj.path
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

think(x::WriteInput, input::Input) = Thunk(x, input)

struct RunCmd{T} <: Action{T}
    calculation::T
end

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
