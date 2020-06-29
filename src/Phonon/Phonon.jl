"""
# module Phonon



# Examples

```jldoctest
julia>
```
"""
module Phonon

using AbInitioSoftwareBase.Inputs: Input, inputstring, write_input

using ..Express: SelfConsistentField, DfptMethod, ForceConstant

export DfptMethod, relay, prep_input, preprocess

function preset end

function relay end

function prep_input end

function preprocess(::DfptMethod, input, template::Input, args...; dry_run = false)
    write_input(input, prep_input(DfptMethod(), template, args...), dry_run)
end

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
