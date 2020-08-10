"""
# module Phonon



# Examples

```jldoctest
julia>
```
"""
module Phonon

using AbInitioSoftwareBase: loadfile
using AbInitioSoftwareBase.Inputs: Input, inputstring, writeinput
using QuantumESPRESSO.Inputs.PWscf: PWInput
using QuantumESPRESSO.Inputs.PHonon: MatdynInput
using QuantumESPRESSO.Outputs.PHonon: parse_frequency, parse_dos
using QuantumESPRESSO.CLI: PWCmd, PhCmd

using ..Express: ElectronicStructure, VibrationalProperty
using ..EosFitting: _check_software_settings

import AbInitioSoftwareBase.Inputs: setcell
import ..Jobs: launchjob

export SelfConsistentField,
    DfptMethod,
    ForceConstant,
    PhononDispersion,
    PhononDensityOfStates,
    prepare,
    launchjob,
    finish,
    load_settings,
    inputstring

struct SelfConsistentField <: ElectronicStructure end
struct DfptMethod <: VibrationalProperty end
struct ForceConstant <: VibrationalProperty end
struct PhononDispersion <: VibrationalProperty end
struct PhononDensityOfStates <: VibrationalProperty end

function setcell(output, template::Input)
    cell = open(output, "r") do io
        str = read(io, String)
        parsecell(str)
    end
    if any(x === nothing for x in cell)
        return  # Not work for relax only
    else
        return setcell(template, cell...)
    end
end

function prepare(::SelfConsistentField, files, outputs, templates; dry_run = false)
    objects = map(files, templates, outputs) do file, template, output
        object = preset_template(SelfConsistentField(), template)
        object = setcell(output, object)
        writeinput(file, object, dry_run)
        object
    end
    return objects
end
prepare(::SelfConsistentField, files, outputs, template::Input; kwargs...) =
    prepare(SelfConsistentField(), files, outputs, fill(template, size(files)))
function prepare(::SelfConsistentField, configfile; kwargs...)
    settings = load_settings(configfile)
    inputs = settings.dirs .* "/phscf.in"
    outputs = settings.dirs .* "/vc-relax.out"
    return prepare(SelfConsistentField(), inputs, outputs, settings.template[1]; kwargs...)
end
function prepare(::DfptMethod, files, templates, args...; dry_run = false)
    objects = map(files, templates, args[1]) do file, template, pw
        object = preset_template(DfptMethod(), template, pw)
        writeinput(file, object, dry_run)
        object
    end
    return objects
end
prepare(::DfptMethod, files, template::Input, args...; kwargs...) =
    prepare(DfptMethod(), files, fill(template, size(files)), args...; kwargs...)
function prepare(calc::DfptMethod, configfile; kwargs...)
    settings = load_settings(configfile)
    inputs = settings.dirs .* "/ph.in"
    pws = map(settings.dirs .* "/phscf.in") do f
        parse(PWInput, read(f, String))
    end
    return prepare(calc, inputs, settings.template[2], pws; kwargs...)
end
function prepare(::ForceConstant, files, templates, args...; dry_run = false)
    objects = map(files, templates) do file, template
        object = preset_template(ForceConstant(), template, args...)
        writeinput(file, object, dry_run)
        object
    end
    return objects
end
prepare(::ForceConstant, files, template::Input, args...; kwargs...) =
    prepare(ForceConstant(), files, fill(template, size(files)), args...; kwargs...)
function prepare(calc::ForceConstant, configfile; kwargs...)
    settings = load_settings(configfile)
    inputs = settings.dirs .* "/q2r.in"
    return prepare(calc, inputs, settings.template[3], settings.template[2]; kwargs...)
end
function prepare(::PhononDispersion, files, templates, args...; dry_run = false)
    objects = map(files, templates) do file, template
        object = preset_template(PhononDispersion(), template, args...)
        writeinput(file, object, dry_run)
        object
    end
    return objects
end
prepare(::PhononDispersion, files, template::Input, args...; kwargs...) =
    prepare(PhononDispersion(), files, fill(template, size(files)), args...; kwargs...)
function prepare(calc::PhononDispersion, configfile; kwargs...)
    settings = load_settings(configfile)
    inputs = settings.dirs .* "/disp.in"
    return prepare(calc, inputs, settings.template[4:-1:2]...; kwargs...)
end

function launchjob(
    ::T,
    configfile;
    kwargs...,
) where {T<:Union{SelfConsistentField,DfptMethod}}
    settings = load_settings(configfile)
    inputs = @. settings.dirs * '/' * (T <: SelfConsistentField ? "scf" : "ph") * ".in"
    outputs = map(Base.Fix2(replace, ".in" => ".out"), inputs)
    return launchjob(
        outputs,
        inputs,
        settings.manager.np,
        settings.bin[T <: SelfConsistentField ? 1 : 2],
        T <: SelfConsistentField ? PWCmd(settings.bin[1]) : PhCmd(settings.bin[2]);
        kwargs...,
    )
end

# function (::Step{PhononDispersion,Prepare{:input}})(
#     matdyn_inputs,
#     q2r_inputs,
#     template::MatdynInput,
# )
#     map(matdyn_inputs, q2r_inputs) do matdyn_input, q2r_input
#         object = parse(Q2rInput, read(InputFile(q2r_input)))
#         template = relay(object, template)
#         if isfile(template.input.flfrq)
#             @set! template.input.flfrq *= if template.input.dos == true  # Phonon DOS calculation
#                 "_dos"  # Append extension `"_dos` to `template.input.flfrq`
#             else  # Phonon dispersion-relation calculation
#                 "_disp"  # Append extension `"_disp` to `template.input.flfrq`
#             end
#         end
#         write(InputFile(matdyn_input), template)
#     end
#     return
# end
# function (::Step{PhononDispersion,Prepare{:input}})(
#     dynmat_inputs,
#     phonon_inputs,
#     template::DynmatInput,
# )
#     map(dynmat_inputs, phonon_inputs) do dynmat_input, phonon_input
#         object = parse(PhInput, read(InputFile(phonon_input)))
#         write(InputFile(dynmat_input), relay(object, template))
#     end
#     return
# end

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

function _expand_settings end

function preset_template end

function parsecell end

function _check_settings(settings)
    map(("template", "pressures", "dir")) do key
        @assert haskey(settings, key)
    end
    @assert isdir(settings["dir"])
    @assert all(isfile.(settings["template"]))
end # function _check_settings

function load_settings(configfile)
    settings = loadfile(configfile)
    _check_settings(settings)  # Errors will be thrown if exist
    return _expand_settings(settings)
end # function load_settings

include("QuantumESPRESSO.jl")

end
