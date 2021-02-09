struct MakeInput{T} <: Action{T} end
function (x::MakeInput{T})(
    file::Union{AbstractString},
    template::Input,
    args...;
    kwargs...,
) where {T<:Union{Scf,LatticeDynamics}}
    input = x(template, args...)
    mkpath(dirname(file))  # In case its parent directory is not created
    writetxt(file, input)
    return input
end
function (x::MakeInput{T})(cfgfile; kwargs...) where {T<:LatticeDynamics}
    settings = loadconfig(cfgfile)
    files = map(dir -> joinpath(dir, shortname(T) * ".in"), settings.dirs)
    previnputs = map(settings.dirs) do dir
        file = joinpath(dir, shortname(prevcalc(T)) * ".in")
        open(file, "r") do io
            parse(inputtype(prevcalc(T)), read(file, String))
        end
    end
    return broadcast(x, files, settings.templates[order(T)], previnputs; kwargs...)
end
function (x::MakeInput{T})(cfgfile; kwargs...) where {T<:Scf}
    settings = loadconfig(cfgfile)
    files = map(dir -> joinpath(dir, shortname(T) * ".in"), settings.dirs)
    templates = settings.templates[1]
    templates = if length(templates) == 1
        fill(templates[1], length(settings.pressures))
    elseif length(templates) != length(settings.pressures)
        throw(DimensionMismatch("!!!"))
    else
        templates
    end
    cells = map(settings.dirs, templates) do dir, template
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
    return broadcast(x, files, templates, first.(cells), last.(cells); kwargs...)
end
function (x::MakeInput{T})(cfgfile; kwargs...) where {T<:Union{PhononDispersion,VDos}}
    settings = loadconfig(cfgfile)
    files = map(dir -> joinpath(dir, shortname(T) * ".in"), settings.dirs)
    ifcinputs = map(settings.dirs) do dir
        file = joinpath(dir, shortname(RealSpaceForceConstants()) * ".in")
        open(file, "r") do io
            parse(inputtype(RealSpaceForceConstants), read(file, String))
        end
    end
    dfptinputs = map(settings.dirs) do dir
        file = joinpath(dir, shortname(Dfpt()) * ".in")
        open(file, "r") do io
            parse(inputtype(Dfpt), read(file, String))
        end
    end
    return broadcast(
        x,
        files,
        settings.templates[order(T)],
        ifcinputs,
        dfptinputs;
        kwargs...,
    )
end

function parsecell end

function inputtype end

buildjob(x::MakeInput, cfgfile) = InternalAtomicJob(() -> x(cfgfile))
