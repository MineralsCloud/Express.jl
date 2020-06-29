module QuantumESPRESSO

using AbInitioSoftwareBase.Inputs: inputstring
using QuantumESPRESSO.Inputs: InputFile, inputstring
using QuantumESPRESSO.Inputs.PWscf:
    AtomicPositionsCard, CellParametersCard, PWInput, optconvert, set_verbosity
using QuantumESPRESSO.Inputs.PHonon: PhInput, Q2rInput, MatdynInput, DynmatInput
using QuantumESPRESSO.Outputs.PWscf: parsefinal
using Setfield: @set!

using ...Express: DfptMethod, SelfConsistentField
import ..Phonon: preset, relay, prep_input, preprocess

function prep_input(::DfptMethod, template::PhInput, from::PWInput)
    template = preset(template)
    return relay(from, template)
end

preprocess(::SelfConsistentField, input, template::PWInput) = write_input(input, template)

# This is a helper function and should not be exported.
function preset(template::PWInput)
    @set! template.control.calculation = "scf"
    return set_verbosity(template, "high")
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

end
