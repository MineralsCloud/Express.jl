"""
# module Phonon



# Examples

```jldoctest
julia>
```
"""
module Phonon

using AbInitioSoftwareBase: load
using AbInitioSoftwareBase.Inputs: Input, writeinput
using SimpleWorkflow: InternalAtomicJob, chain

import ..Express
using ..Express:
    Calculation,
    LatticeDynamics,
    SelfConsistentField,
    Scf,
    FixedIonSelfConsistentField,
    Action,
    MakeCmd,
    distprocs,
    makescript,
    load_settings
using ..EosFitting: VcOptim

import AbInitioSoftwareBase.Inputs: set_cell

export Dfpt,
    Dfpt,
    SelfConsistentField,
    Scf,
    FixedIonSelfConsistentField,
    RealSpaceForceConstants,
    PhononDispersion,
    PhononDensityOfStates,
    VDos,
    ZoneCenterPhonons,
    ZoneCentrePhonons,
    MakeInput,
    makescript,
    writeinput,
    standardize,
    standardize,
    customize,
    load_settings

struct Dfpt <: LatticeDynamics end
struct RealSpaceForceConstants <: LatticeDynamics end
struct PhononDispersion <: LatticeDynamics end
struct PhononDensityOfStates <: LatticeDynamics end
struct ZoneCenterPhonons <: LatticeDynamics end
const VDos = PhononDensityOfStates
const ZoneCentrePhonons = ZoneCenterPhonons  # For English users

function set_cell(template::Input, output)
    cell = open(output, "r") do io
        str = read(io, String)
        parsecell(str)
    end
    if any(x === nothing for x in cell)
        error("set cell failed!")
    else
        return set_cell(template, cell...)
    end
end

struct MakeInput{T} <: Action{T} end
MakeInput(::T) where {T<:Calculation} = MakeInput{T}()
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
    settings = load_settings(cfgfile)
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
    settings = load_settings(cfgfile)
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
    settings = load_settings(cfgfile)
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

buildjob(::MakeInput{T}, cfgfile) where {T} =
    InternalAtomicJob(() -> MakeInput(T())(cfgfile))
function buildjob(x::MakeCmd{T}, cfgfile) where {T}
    settings = load_settings(cfgfile)
    inp = map(dir -> joinpath(dir, shortname(T) * ".in"), settings.dirs)
    out = map(dir -> joinpath(dir, shortname(T) * ".out"), settings.dirs)
    return Express.buildjob(
        x,
        out,
        inp,
        settings.manager.np,
        settings.bin[order(T)];
        use_shell = settings.use_shell,
    )
end

function buildworkflow(cfgfile)
    step1 = buildjob(MakeInput(Scf()), cfgfile)
    step2 = chain(step1, buildjob(MakeCmd(Scf()), cfgfile)[1])
    step3 = chain(step2[end], buildjob(MakeInput(Dfpt()), cfgfile))
    step4 = chain(step3[end], buildjob(MakeCmd(Dfpt()), cfgfile)[1])
    step5 = chain(step4[end], buildjob(MakeInput(RealSpaceForceConstants()), cfgfile))
    step6 = chain(step5[end], buildjob(MakeCmd(RealSpaceForceConstants()), cfgfile)[1])
    settings = load(cfgfile)
    x = if settings["workflow"] == "phonon dispersion"
        PhononDispersion()
    elseif settings["workflow"] == "vdos"
        VDos()
    else
        error("unsupported option!")
    end
    step7 = chain(step6[end], buildjob(MakeInput(x), cfgfile))
    step8 = chain(step7[end], buildjob(MakeCmd(x), cfgfile)[1])
    return step8
end

order(x) = order(typeof(x))
order(::Type{Scf}) = 1
order(::Type{Dfpt}) = 2
order(::Type{RealSpaceForceConstants}) = 3
order(::Type{PhononDispersion}) = 4
order(::Type{PhononDensityOfStates}) = 4

prevcalc(x) = prevcalc(typeof(x))
prevcalc(::Type{Scf}) = VcOptim()
prevcalc(::Type{Dfpt}) = Scf()
prevcalc(::Type{RealSpaceForceConstants}) = Dfpt()
prevcalc(::Type{PhononDispersion}) = RealSpaceForceConstants()
prevcalc(::Type{PhononDensityOfStates}) = RealSpaceForceConstants()

shortname(x::Calculation) = shortname(typeof(x))

function standardize end

function customize end

function parseoutput end

function expand_settings end

function check_software_settings end

function inputtype end

function parsecell end

function check_settings(settings)
    map(("templates", "pressures", "workdir")) do key
        @assert haskey(settings, key)
    end
    if !isdir(expanduser(settings["workdir"]))
        @warn "`workdir` is not reachable, be careful!"
    end
    for paths in settings["templates"]
        for path in paths
            if !isfile(path)
                @warn "template \"$path\" is not reachable, be careful!"
            end
        end
    end
end

end
