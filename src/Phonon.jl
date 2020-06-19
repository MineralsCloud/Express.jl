"""
# module Phonon



# Examples

```jldoctest
julia>
```
"""
module Phonon

using ConstructionBase: setproperties
using QuantumESPRESSO.Inputs: InputFile, qestring
using QuantumESPRESSO.Inputs.PWscf:
    AtomicPositionsCard, CellParametersCard, PWInput, optconvert
using QuantumESPRESSO.Inputs.PHonon: PhInput, Q2rInput, MatdynInput, DynmatInput
using QuantumESPRESSO.Outputs: OutputFile
using QuantumESPRESSO.Outputs.PWscf: parsefinal
using QuantumESPRESSO.Setters: CellParametersSetter, VerbositySetter
using Setfield: @set!

using ..Express:
    Calculation,
    SelfConsistentField,
    DfptMethod,
    ForceConstant,
    PREPARE_INPUT,
    LAUNCH_JOB,
    ANALYSE_OUTPUT

export set_structure, relay, preprocess

"""
    set_structure(output::AbstractString, template::PWInput)

Read structure information from `output`, and update related fields of `template`.
"""
function set_structure(output::AbstractString, template::PWInput)
    str = read(OutputFile(output))
    cell = parsefinal(CellParametersCard{Float64}, str)
    return setproperties(
        template,
        atomic_positions = parsefinal(AtomicPositionsCard, str),
        cell_parameters = CellParametersCard(cell.data / template.system.celldm[1], "alat"), # The result of `parsefinal` must be a `CellParametersCard` with `"bohr"` or `"angstrom"` option, convert it to "bohr" by default
    )
end # function set_structure

# This is a helper function and should not be exported.
function preset(template::PWInput)
    @set! template.control.verbosity = "high"
    @set! template.control.wf_collect = true
    @set! template.control.tstress = true
    @set! template.control.tprnfor = true
    @set! template.control.disk_io = "high"
    @set! template.control.calculation = "scf"
    return template
end # function preset
function preset(template::PhInput)
    @set! template.inputph.verbosity = "high"
    return template
end # function preset

"""
    relay(from::PWInput, to::PhInput)

Relay shared information from a `PWInput` to a `PhInput`.

A `PWInput` before a `PhInput` has the information of `outdir` and `prefix`. They must keep the same in a
phonon calculation.
"""
function relay(pw::PWInput, ph::PhInput)
    @set! ph.inputph.outdir = pw.control.outdir
    @set! ph.inputph.prefix = pw.control.prefix
    return ph
end # function relay
"""
    relay(from::PhInput, to::Q2rInput)

Relay shared information from a `PhInput` to a `Q2rInput`.

A `PhInput` before a `Q2rInput` has the information of `fildyn`. It must keep the same in a q2r calculation.
"""
function relay(ph::PhInput, q2r::Q2rInput)
    @set! q2r.input.fildyn = ph.inputph.fildyn
    return q2r
end # function relay
"""
    relay(from::Q2rInput, to::MatdynInput)

Relay shared information from a `Q2rInput` to a `MatdynInput`.

A `Q2rInput` before a `MatdynInput` has the information of `fildyn`, `flfrc` and `loto_2d`. They must keep the same
in a matdyn calculation.
"""
function relay(q2r::Q2rInput, matdyn::MatdynInput)
    @set! matdyn.input.flfrc = q2r.input.flfrc
    @set! matdyn.input.loto_2d = q2r.input.loto_2d
    return matdyn
end # function relay
"""
    relay(from::PhInput, to::DynmatInput)

Relay shared information from a `PhInput` to a `DynmatInput`.

A `PhInput` before a `DynmatInput` has the information of `asr`, `fildyn` and `amass`. They must keep the same
in a dynmat calculation.
"""
function relay(ph::PhInput, dynmat::DynmatInput)
    # @set! dynmat.input.asr = ph.inputph.asr  # TODO
    @set! dynmat.input.fildyn = ph.inputph.fildyn
    @set! dynmat.input.amass = ph.inputph.amass
    return dynmat
end # function relay

function (::Step{SelfConsistentField,PREPARE_INPUT})(
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
function preprocess(
    ::Step{DfptMethod,PREPARE_INPUT},
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
function preprocess(
    ::Step{ForceConstant,PREPARE_INPUT},
    q2r_inputs,
    phonon_inputs,
    template::Q2rInput,
)
    map(q2r_inputs, phonon_inputs) do q2r_input, phonon_input
        object = parse(PhInput, read(InputFile(phonon_input)))
        write(InputFile(q2r_input), relay(object, template))
    end
    return
end # function preprocess
function preprocess(
    ::Step{PhononDispersion,PREPARE_INPUT},
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
function preprocess(
    ::Step{PhononDispersion,PREPARE_INPUT},
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

#???

end
