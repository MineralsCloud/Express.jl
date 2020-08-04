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

using ..Express: SelfConsistentField, DfptMethod, ForceConstant, PhononDispersion
using ..EosFitting: _check_software_settings

import AbInitioSoftwareBase.Inputs: set_structure
import ..Express

export SelfConsistentField,
    DfptMethod,
    ForceConstant,
    PhononDispersion,
    prepare,
    process,
    load_settings,
    inputstring

function set_structure(output, template::Input)
    cell = open(output, "r") do io
        str = read(io, String)
        parsecell(str)
    end
    if any(x === nothing for x in cell)
        return  # Not work for relax only
    else
        return set_structure(template, cell...)
    end
end

function prepare(::SelfConsistentField, files, outputs, templates; dry_run = false)
    objects = map(files, templates, outputs) do file, template, output
        object = preset_template(SelfConsistentField(), template)
        object = set_structure(output, object)
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
    objects = map(files, templates) do file, template
        object = preset_template(DfptMethod(), template, args...)
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
    return prepare(calc, inputs, settings.template[2], settings.template[1]; kwargs...)
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
        T isa SelfConsistentField ? PWCmd() : PhCmd();
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

function load_settings(path)
    settings = loadfile(path)
    _check_settings(settings)  # Errors will be thrown if exist
    return _expand_settings(settings)
end # function load_settings

include("QuantumESPRESSO.jl")

end
