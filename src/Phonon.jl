"""
# module Phonon



# Examples

```jldoctest
julia>
```
"""
module Phonon

using AbInitioSoftwareBase.CLI: Mpiexec
using AbInitioSoftwareBase.Inputs: Input, writeinput
using SimpleWorkflow: ExternalAtomicJob, InternalAtomicJob, chain, parallel

using ..Express:
    Calculation,
    LatticeDynamics,
    SelfConsistentField,
    Scf,
    FixedIonSelfConsistentField,
    distprocs,
    makescript,
    load_settings
using ..EosFitting: VcOptim

import AbInitioSoftwareBase.Inputs: set_cell

export DensityFunctionalPerturbationTheory,
    Dfpt,
    SelfConsistentField,
    Scf,
    FixedIonSelfConsistentField,
    InteratomicForceConstants,
    Ifc,
    PhononDispersion,
    PhononDensityOfStates,
    VDos,
    ZoneCenterPhonons,
    ZoneCentrePhonons,
    standardize,
    makeinput,
    makescript,
    writeinput,
    standardize,
    customize,
    load_settings,
    buildjob

struct DensityFunctionalPerturbationTheory <: LatticeDynamics end
struct InteratomicForceConstants <: LatticeDynamics end
struct PhononDispersion <: LatticeDynamics end
struct PhononDensityOfStates <: LatticeDynamics end
struct ZoneCenterPhonons <: LatticeDynamics end
const Dfpt = DensityFunctionalPerturbationTheory
const Ifc = InteratomicForceConstants
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

struct MakeInput{T<:Calculation}
    calc::T
end
function (x::MakeInput{<:Union{Scf,LatticeDynamics}})(
    file,
    template::Input,
    args...;
    kwargs...,
)
    object = customize(standardize(template, x.calc), args...; kwargs...)
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
function (x::MakeInput{<:LatticeDynamics})(cfgfile; kwargs...)
    calc = x.calc
    settings = load_settings(cfgfile)
    files = map(dir -> joinpath(dir, shortname(calc) * ".in"), settings.dirs)
    previnputs = map(settings.dirs) do dir
        file = joinpath(dir, shortname(prevcalc(calc)) * ".in")
        open(file, "r") do io
            parse(inputtype(prevcalc(calc)), read(file, String))
        end
    end
    return x(files, settings.templates[order(calc)], previnputs; kwargs...)
end
function (x::MakeInput{Scf})(cfgfile; kwargs...)
    calc = x.calc
    settings = load_settings(cfgfile)
    files = map(dir -> joinpath(dir, shortname(calc) * ".in"), settings.dirs)
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
function (x::MakeInput{<:Union{PhononDispersion,VDos}})(cfgfile; kwargs...)
    calc = x.calc
    settings = load_settings(cfgfile)
    files = map(dir -> joinpath(dir, shortname(calc) * ".in"), settings.dirs)
    ifcinputs = map(settings.dirs) do dir
        file = joinpath(dir, shortname(Ifc()) * ".in")
        open(file, "r") do io
            parse(inputtype(Ifc()), read(file, String))
        end
    end
    dfptinputs = map(settings.dirs) do dir
        file = joinpath(dir, shortname(Dfpt()) * ".in")
        open(file, "r") do io
            parse(inputtype(Dfpt()), read(file, String))
        end
    end
    return x(files, settings.templates[order(calc)], ifcinputs, dfptinputs; kwargs...)
end

makeinput(calc::Calculation) = MakeInput(calc)

function buildjob(x::MakeInput)
    function _buildjob(cfgfile)
        InternalAtomicJob(() -> x(cfgfile))
    end
end
function buildjob(calc::Calculation)
    function _buildjob(outputs, inputs, np, exe; kwargs...)
        # `map` guarantees they are of the same size, no need to check.
        n = distprocs(np, length(inputs))
        subjobs = map(outputs, inputs) do output, input
            f = Mpiexec(n; kwargs...) ∘ exe
            cmd = f(stdin = input, stdout = output)
            ExternalAtomicJob(cmd)
        end
        return parallel(subjobs...)
    end
    function _buildjob(template, view)
        ExternalAtomicJob(makescript(template, view))
    end
    function _buildjob(cfgfile)
        settings = load_settings(cfgfile)
        inp = map(dir -> joinpath(dir, shortname(calc) * ".in"), settings.dirs)
        out = map(dir -> joinpath(dir, shortname(calc) * ".out"), settings.dirs)
        return _buildjob(out, inp, settings.manager.np, settings.bin[order(calc)])
    end
end

function buildworkflow(cfgfile)
    step1 = buildjob(makeinput(Scf()))(cfgfile)
    step12 = chain(step1, buildjob(Scf())(cfgfile)[1])
    step123 = chain(step12[end], buildjob(makeinput(Dfpt()))(cfgfile))
    step1234 = chain(step123[end], buildjob(Dfpt())(cfgfile)[1])
    step12345 = chain(step1234[end], buildjob(makeinput(Ifc()))(cfgfile))
    step123456 = chain(step12345[end], buildjob(Ifc())(cfgfile)[1])
    step1234567 = chain(step123456[end], buildjob(makeinput(PhononDispersion()))(cfgfile))
    step12345678 = chain(step1234567[end], buildjob(PhononDispersion())(cfgfile)[1])
    return step12345678
end

order(::Scf) = 1
order(::Dfpt) = 2
order(::Ifc) = 3
order(::PhononDispersion) = 4
order(::PhononDensityOfStates) = 4

prevcalc(::Scf) = VcOptim()
prevcalc(::Dfpt) = Scf()
prevcalc(::Ifc) = Dfpt()
prevcalc(::PhononDispersion) = Ifc()
prevcalc(::PhononDensityOfStates) = Ifc()

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
