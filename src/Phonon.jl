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
using QuantumESPRESSO.Inputs.PHonon: MatdynInput
using QuantumESPRESSO.Outputs.PHonon: parse_frequency, parse_dos
using QuantumESPRESSO.CLI: PWCmd, PhCmd

using ..Express: ElectronicStructure, VibrationalProperty, Scf

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

function makeinput(calc::Union{Scf,Dfpt,Ifc,PhononDispersion,VDos})
    function _makeinput(file, template::Input, args...; kwargs...)
        object = customize(standardize(template, calc), args...; kwargs...)
        writeinput(file, object)
        return object
    end
    function _makeinput(files, templates, restargs...; kwargs...)
        if templates isa Input
            templates = fill(templates, size(files))
        end
        objects = map(files, templates, zip(restargs...)) do file, template, args
            _makeinput(file, template, args...; kwargs...)
        end
        return objects
    end
    function _makeinput(cfgfile; kwargs...)
        settings = load_settings(cfgfile)
        files = map(dir -> joinpath(dir, shortname(calc) * ".in"), settings.dirs)
        previnputs = map(settings.dirs) do dir
            file = joinpath(dir, shortname(calc) * ".in")
            open(file, "r") do io
                parse(previnputtype(calc), read(file, String))
            end
        end
        return _makeinput(files, settings.templates, previnputs; kwargs...)
    end
end

function launchjob(::T, configfile; kwargs...) where {T<:Union{Scf,Dfpt}}
    settings = load_settings(configfile)
    inputs = @. settings.dirs * '/' * (T <: Scf ? "scf" : "ph") * ".in"
    outputs = map(Base.Fix2(replace, ".in" => ".out"), inputs)
    return launchjob(
        outputs,
        inputs,
        settings.manager.np,
        settings.bin[T <: Scf ? 1 : 2],
        T <: Scf ? PWCmd(settings.bin[1]) : PhCmd(settings.bin[2]);
        kwargs...,
    )
end

function finish(::PhononDispersion, outputs::AbstractArray)
    map(outputs) do output
        str = read(output, String)
        parse_frequency(str)
    end
end
function finish(::PhononDispersion, configfile)
    settings = load_settings(configfile)
    inputs = settings.dirs .* "/disp.in"
    outputs = map(inputs, settings.dirs) do input, dir
        joinpath(dir, parse(MatdynInput, read(input, String)).input.flfrq)
    end
    return finish(PhononDispersion(), outputs)
end
function finish(::PhononDensityOfStates, outputs::AbstractArray)
    map(outputs) do output
        str = read(output, String)
        parse_dos(str)
    end
end
function finish(::PhononDensityOfStates, configfile)
    settings = load_settings(configfile)
    inputs = settings.dirs .* "/disp.in"
    outputs = map(inputs, settings.dirs) do input, dir
        joinpath(dir, parse(MatdynInput, read(input, String)).input.fldos)
    end
    return finish(PhononDensityOfStates(), outputs)
end

function standardize end

function customize end

function parseoutput end

function expand_settings end

function check_software_settings end

function shortname end

function previnputtype end

function preset_template end

function parsecell end

function _check_settings(settings)
    map(("template", "pressures", "dir")) do key
        @assert haskey(settings, key)
    end
    @assert isdir(settings["dir"])
    @assert all(isfile.(settings["template"]))
end

function load_settings(configfile)
    settings = load(configfile)
    _check_settings(settings)  # Errors will be thrown if exist
    return expand_settings(settings)
end

end
