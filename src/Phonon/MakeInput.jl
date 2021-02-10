struct MakeInput{T} <: Action{T} end
function (x::MakeInput{T})(
    file::Union{AbstractString},
    template::Input,
    args...,
) where {T<:Union{Scf,LatticeDynamics}}
    input = x(template, args...)
    mkpath(dirname(file))  # In case its parent directory is not created
    writetxt(file, input)
    return input
end
function (x::MakeInput{T})(cfgfile) where {T<:LatticeDynamics}
    config = loadconfig(cfgfile)
    files = map(dir -> joinpath(dir, shortname(T) * ".in"), config.dirs)
    previnputs = map(config.dirs) do dir
        file = joinpath(dir, shortname(prevcalc(T)) * ".in")
        open(file, "r") do io
            parse(inputtype(prevcalc(T)), read(file, String))
        end
    end
    return broadcast(
        x,
        files,
        map(o -> getfield(o, order(T)), config.templates),
        previnputs,
    )
end
function (x::MakeInput{T})(cfgfile) where {T<:Scf}
    config = loadconfig(cfgfile)
    files = map(dir -> joinpath(dir, shortname(T) * ".in"), config.dirs)
    cells = map(config.dirs) do dir
        file = joinpath(dir, shortname(VcOptim()) * ".out")
        cell = open(file, "r") do io
            str = read(io, String)
            parsecell(str)
        end
        if any(x === nothing for x in cell)
            error("set cell failed!")
        else
            cell
        end
    end
    return broadcast(
        x,
        files,
        map(o -> getfield(o, order(T)), config.templates),
        first.(cells),
        last.(cells),
    )
end
function (x::MakeInput{T})(cfgfile) where {T<:Union{PhononDispersion,VDos}}
    config = loadconfig(cfgfile)
    files = map(dir -> joinpath(dir, shortname(T) * ".in"), config.dirs)
    ifcinputs = map(config.dirs) do dir
        file = joinpath(dir, shortname(RealSpaceForceConstants()) * ".in")
        open(file, "r") do io
            parse(inputtype(RealSpaceForceConstants), read(file, String))
        end
    end
    dfptinputs = map(config.dirs) do dir
        file = joinpath(dir, shortname(Dfpt()) * ".in")
        open(file, "r") do io
            parse(inputtype(Dfpt), read(file, String))
        end
    end
    return broadcast(
        x,
        files,
        map(o -> getfield(o, order(T)), config.templates),
        ifcinputs,
        dfptinputs,
    )
end

function parsecell end

function inputtype end

buildjob(x::MakeInput, cfgfile) = InternalAtomicJob(() -> x(cfgfile))
