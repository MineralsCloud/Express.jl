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
    return broadcast(x, files, config.templates[order(T)], previnputs)
end
function (x::MakeInput{T})(cfgfile) where {T<:Scf}
    config = loadconfig(cfgfile)
    files = map(dir -> joinpath(dir, shortname(T) * ".in"), config.dirs)
    templates = config.templates[1]
    templates = if length(templates) == 1
        fill(templates[1], length(config.pressures))
    elseif length(templates) != length(config.pressures)
        throw(DimensionMismatch("!!!"))
    else
        templates
    end
    cells = map(config.dirs, templates) do dir, template
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
    return broadcast(x, files, templates, first.(cells), last.(cells))
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
    return broadcast(x, files, config.templates[order(T)], ifcinputs, dfptinputs)
end

function parsecell end

function inputtype end

buildjob(x::MakeInput, cfgfile) = InternalAtomicJob(() -> x(cfgfile))
