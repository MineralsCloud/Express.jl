"""
# module Phonon



# Examples

```jldoctest
julia>
```
"""
module Phonon

using AbInitioSoftwareBase: load
using AbInitioSoftwareBase.Inputs: Input, inputstring, writeinput

using ..Express: Calculation, ElectronicStructure, VibrationalProperty, Scf
using ..EosFitting: VariableCellOptimization

import AbInitioSoftwareBase.Inputs: set_cell

struct DensityFunctionalPerturbationTheory <: VibrationalProperty end
struct InteratomicForceConstants <: VibrationalProperty end
struct PhononDispersion <: VibrationalProperty end
struct PhononDensityOfStates <: VibrationalProperty end
const Dfpt = DensityFunctionalPerturbationTheory
const Ifc = InteratomicForceConstants
const VDos = PhononDensityOfStates

function set_cell(output, template::Input)
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
function (x::MakeInput{<:VibrationalProperty})(file, template::Input, args...; kwargs...)
    object = customize(standardize(template, x.calc), args...; kwargs...)
    writeinput(file, object)
    return object
end
function (x::MakeInput{<:VibrationalProperty})(files, templates, restargs...; kwargs...)
    if templates isa Input
        templates = fill(templates, size(files))
    end
    objects = map(files, templates, zip(restargs...)) do file, template, args
        x(file, template, args...; kwargs...)
    end
    return objects
end
function (x::MakeInput{<:VibrationalProperty})(cfgfile; kwargs...)
    calc = x.calc
    settings = load_settings(cfgfile)
    files = map(dir -> joinpath(dir, shortname(calc) * ".in"), settings.dirs)
    previnputs = map(settings.dirs) do dir
        file = joinpath(dir, shortname(calc) * ".in")
        open(file, "r") do io
            parse(previnputtype(calc), read(file, String))
        end
    end
    return x(files, settings.templates, previnputs...; kwargs...)
end
function (x::MakeInput{Scf})(cfgfile; kwargs...)
    calc = x.calc
    settings = load_settings(cfgfile)
    files = map(dir -> joinpath(dir, shortname(calc) * ".in"), settings.dirs)
    newcells = map(settings.dirs) do dir
        file = joinpath(dir, shortname(VariableCellOptimization()) * ".out")
        open(file, "r") do io
            parsecell(read(file, String))
        end
    end
    return x(files, settings.templates, newcells; kwargs...)
end

makeinput(calc::Calculation) = MakeInput(calc)

function standardize end

function customize end

function parseoutput end

function expand_settings end

function check_software_settings end

function shortname end

function previnputtype end

function parsecell end

function check_settings(settings)
    map(("templates", "pressures", "workdir")) do key
        @assert haskey(settings, key)
    end
    @assert isdir(settings["workdir"])
    @assert all(map(isfile, settings["templates"]))
end

function load_settings(configfile)
    settings = load(configfile)
    check_settings(settings)  # Errors will be thrown if exist
    return expand_settings(settings)
end

end
