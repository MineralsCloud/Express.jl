"""
# module Phonon



# Examples

```jldoctest
julia>
```
"""
module Phonon

using AbInitioSoftwareBase.Inputs: inputstring

using ..Express: Step, Calculation, SelfConsistentField, DfptMethod, ForceConstant, Prepare

export relay, Step

function preset end

function relay end

function (::Step{DfptMethod,Prepare{:input}})(template::PhInput, from::PWInput)
    template = preset(template)
    return relay(from, template)
end
function (::Step{DfptMethod,Prepare{:input}})(
    input,
    template::PhInput,
    from::PWInput;
    dry_run = false,
)
    template = preset(template)
    object = relay(from, template)
    if dry_run
        if isfile(input)
            @warn "file `$input` will be overwritten!"
        else
            @warn "file `$input` will be created!"
        end
        print(inputstring(object))
    else
        mkpath(dirname(input))
        open(input, "w") do io
            write(io, inputstring(object))
        end
    end
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

end
