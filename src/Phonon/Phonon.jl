"""
# module Phonon



# Examples

```jldoctest
julia>
```
"""
module Phonon

using ..Express: Step, Calculation, SelfConsistentField, DfptMethod, ForceConstant, Prepare

export relay, Step

function preset end

function relay end

function (::Step{SelfConsistentField,Prepare{:input}})(
    scf_inputs,
    vc_outputs,
    template::PWInput,
)
    template = preset(template)
    map(scf_inputs, vc_outputs) do input, output
        # Get a new `object` from the `template`, with its `alat` and `pressure` changed
        object = set_structure(output, template)
        write(InputFile(input), object)
        object
    end
end
function (::Step{DfptMethod,Prepare{:input}})(
    phonon_inputs,
    pwscf_inputs,
    template::PhInput,
)
    template = preset(template)
    map(phonon_inputs, pwscf_inputs) do phonon_input, pwscf_input
        object = parse(PWInput, read(InputFile(pwscf_input)))
        write(InputFile(phonon_input), relay(object, template))
    end
    return
end
function (::Step{ForceConstant,Prepare{:input}})(
    q2r_inputs,
    phonon_inputs,
    template::Q2rInput,
)
    map(q2r_inputs, phonon_inputs) do q2r_input, phonon_input
        object = parse(PhInput, read(InputFile(phonon_input)))
        write(InputFile(q2r_input), relay(object, template))
    end
    return
end
function (::Step{PhononDispersion,Prepare{:input}})(
    matdyn_inputs,
    q2r_inputs,
    template::MatdynInput,
)
    map(matdyn_inputs, q2r_inputs) do matdyn_input, q2r_input
        object = parse(Q2rInput, read(InputFile(q2r_input)))
        template = relay(object, template)
        if isfile(template.input.flfrq)
            @set! template.input.flfrq *= if template.input.dos == true  # Phonon DOS calculation
                "_dos"  # Append extension `"_dos` to `template.input.flfrq`
            else  # Phonon dispersion-relation calculation
                "_disp"  # Append extension `"_disp` to `template.input.flfrq`
            end
        end
        write(InputFile(matdyn_input), template)
    end
    return
end
function (::Step{PhononDispersion,Prepare{:input}})(
    dynmat_inputs,
    phonon_inputs,
    template::DynmatInput,
)
    map(dynmat_inputs, phonon_inputs) do dynmat_input, phonon_input
        object = parse(PhInput, read(InputFile(phonon_input)))
        write(InputFile(dynmat_input), relay(object, template))
    end
    return
end

end
