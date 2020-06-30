"""
# module Phonon



# Examples

```jldoctest
julia>
```
"""
module Phonon

using AbInitioSoftwareBase.Inputs: Input, inputstring, write_input
using OptionalArgChecks: @argcheck

using ..Express: SelfConsistentField, DfptMethod, ForceConstant, load_settings
using ..EosFitting: _check_software_settings
import ..Express

export DfptMethod, relay, prep_input, preprocess, load_settings

function preset end

function relay end

function prep_input end

function preprocess(::DfptMethod, inputs, template::Input, args...; dry_run = false)
    map(inputs) do input
        write_input(input, prep_input(DfptMethod(), template, args...), dry_run)
    end
end
function preprocess(calc::DfptMethod, path; kwargs...)
    settings = load_settings(path)
    inputs = settings.dirs .* "/ph.in"
    return preprocess(calc, inputs, settings.template[2], settings.template[1]; kwargs...)
end

function Express._check_settings(settings)
    map(("template", "pressures", "dir")) do key
        @argcheck haskey(settings, key)
    end
    @assert isdir(settings["dir"])
    @assert all(isfile.(settings["template"]))
end # function _check_settings

# function (::Step{ForceConstant,Prepare{:input}})(
#     q2r_inputs,
#     phonon_inputs,
#     template::Q2rInput,
# )
#     map(q2r_inputs, phonon_inputs) do q2r_input, phonon_input
#         object = parse(PhInput, read(InputFile(phonon_input)))
#         write(InputFile(q2r_input), relay(object, template))
#     end
#     return
# end

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

include("QuantumESPRESSO.jl")

end
