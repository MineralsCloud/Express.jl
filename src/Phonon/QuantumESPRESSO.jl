module QuantumESPRESSO

using Crystallography: Cell, eachatom
using QuantumESPRESSO.Inputs: InputFile, inputstring
using QuantumESPRESSO.Inputs.PWscf:
    AtomicPositionsCard, CellParametersCard, PWInput, optconvert
using QuantumESPRESSO.Inputs.PHonon: PhInput, Q2rInput, MatdynInput, DynmatInput
using QuantumESPRESSO.Outputs.PWscf: parsefinal
using Setfield: @set!

import Express.Phonon: set_structure, preset

function set_structure(
    template::PWInput,
    cell_parameters::CellParametersCard,
    atomic_positions::Union{Nothing,AtomicPositionsCard} = nothing,
)
    @set! template.atomic_positions = atomic_positions
    @set! template.cell_parameters = cell_parameters
    return template
end # function set_structure
function set_structure(template::PWInput, cell::Cell)
    @set! template.atomic_positions = AtomicPositionsCard([atom for atom in eachatom(cell)])
    @set! template.cell_parameters = CellParametersCard(cell.lattice)
    return template
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

function (::Step{SelfConsistentField,PREPARE_INPUT})(vc_output, template::PWInput)
    str = open(vc_output, "r") do io
        read(vc_output, String)
    end
    cell = parsefinal(CellParametersCard{Float64}, str)
    cell_parameters =
        CellParametersCard(cell.data / template.system.celldm[1], "alat"), # The result of `parsefinal` must be a `CellParametersCard` with `"bohr"` or `"angstrom"` option, convert it to "bohr" by default
        atomic_positions = tryparsefinal(AtomicPositionsCard, str)
    set_structure(template, cell_parameters, atomic_positions)
    return preset(template)
end

end
