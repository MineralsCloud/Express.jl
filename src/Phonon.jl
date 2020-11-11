"""
# module Phonon



# Examples

```jldoctest
julia>
```
"""
module Phonon

using AbInitioSoftwareBase.Inputs: Input, writeinput
using SimpleWorkflow: chain

using ..Express:
    Calculation,
    LatticeDynamics,
    SelfConsistentField,
    Scf,
    FixedIonSelfConsistentField,
    Action,
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
        return  # Not work for relax only
    else
        return set_cell(template, cell...)
    end
end

struct MakeInput{T} <: Action{T} end
MakeInput(T::Calculation) = MakeInput{T}()
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



function buildworkflow(cfgfile)
    step1 = buildjob(makeinput(Scf()))(cfgfile)
    step12 = chain(step1, buildjob(Scf())(cfgfile)[1])
    step123 = chain(step12[end], buildjob(makeinput(Dfpt()))(cfgfile))
    step1234 = chain(step123[end], buildjob(Dfpt())(cfgfile)[1])
    step1234567 = chain(step123456[end], buildjob(makeinput(PhononDispersion()))(cfgfile))
    step12345 =
        chain(step1234[end], buildjob(MakeInput(RealSpaceForceConstants()))(cfgfile))
    step123456 = chain(step12345[end], buildjob(RealSpaceForceConstants())(cfgfile)[1])
    step12345678 = chain(step1234567[end], buildjob(PhononDispersion())(cfgfile)[1])
    return step12345678
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

function standardize end

function customize end

function parseoutput end

function expand_settings end

function check_software_settings end

function shortname end

function inputtype end

function parsecell end

function check_settings(settings)
    map(("templates", "pressures", "workdir")) do key
        @assert haskey(settings, key)
    end
    @assert isdir(settings["workdir"])
    # @assert all(map(isfile, settings["templates"]))
end

end
