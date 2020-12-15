struct MakeInput{T} <: Action{T} end
function (x::MakeInput{T})(
    file,
    template::Input,
    args...;
    kwargs...,
) where {T<:Union{Scf,LatticeDynamics}}
    object = customize(standardize(template, T()), args...; kwargs...)
    writeinput(file, object)
    return object
end
function (x::MakeInput{<:Union{Scf,LatticeDynamics}})(
    files,
    templates,
    restargs...;
    kwargs...,
)
    if templates isa Input
        templates = fill(templates, size(files))
    end
    objects = map(files, templates, zip(restargs...)) do file, template, args
        x(file, template, args...; kwargs...)
    end
    return objects
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
    return x(files, settings.templates[order(T)], previnputs; kwargs...)
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
    templates = map(settings.dirs, templates) do dir, template
        file = joinpath(dir, shortname(VcOptim()) * ".out")
        set_cell(template, file)
    end
    return x(files, templates; kwargs...)
end
function (x::MakeInput{T})(cfgfile; kwargs...) where {T<:Union{PhononDispersion,VDos}}
    settings = loadconfig(cfgfile)
    files = map(dir -> joinpath(dir, shortname(T) * ".in"), settings.dirs)
    ifcinputs = map(settings.dirs) do dir
        file = joinpath(dir, shortname(RealSpaceForceConstants()) * ".in")
        open(file, "r") do io
            parse(inputtype(RealSpaceForceConstants()), read(file, String))
        end
    end
    dfptinputs = map(settings.dirs) do dir
        file = joinpath(dir, shortname(Dfpt()) * ".in")
        open(file, "r") do io
            parse(inputtype(Dfpt()), read(file, String))
        end
    end
    return x(files, settings.templates[order(T)], ifcinputs, dfptinputs; kwargs...)
end

function standardize end

function customize end

buildjob(::MakeInput{T}, cfgfile) where {T} =
    InternalAtomicJob(() -> MakeInput(T())(cfgfile))
