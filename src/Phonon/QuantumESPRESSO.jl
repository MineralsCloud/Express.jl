module QuantumESPRESSO

using AbInitioSoftwareBase.Inputs: inputstring, writeinput
using Distributed: LocalManager
using QuantumESPRESSO.Inputs.PWscf:
    AtomicPositionsCard, CellParametersCard, PWInput, optconvert, set_verbosity
using QuantumESPRESSO.Inputs.PHonon: PhInput, Q2rInput, MatdynInput, DynmatInput
using QuantumESPRESSO.Outputs.PWscf: parsefinal
using Setfield: @set!
using Unitful: @u_str
using UnitfulAtomic

using ...Express: DfptMethod, SelfConsistentField, ForceConstant
import ...EosFitting: parsecell
import ..Phonon: preset, relay, prep_input, preprocess, _expand_settings, Step

function prep_input(::DfptMethod, template::PhInput, from::PWInput)
    template = preset(template)
    return relay(from, template)
end

prepare(::SelfConsistentField, input, template::PWInput) = writeinput(input, template)

# function (::Step{ForceConstant,Action{:prepare_input}})(
#     q2r_input,
#     phonon_input,
#     template::Q2rInput,
# )
#     object = parse(PhInput, read(phonon_input, String))
#     write(q2r_input, inputstring(relay(object, template)))
#     return
# end

# This is a helper function and should not be exported.
function preset(template::PWInput)
    @set! template.control.calculation = "scf"
    return set_verbosity(template, "high")
end # function preset
function preset(template::PhInput)
    @set! template.inputph.verbosity = "high"
    return template
end # function preset

function _expand_settings(settings)
    templatetexts = [read(expanduser(f), String) for f in settings["template"]]
    template = parse(PWInput, templatetexts[1]), parse(PhInput, templatetexts[2]), parse(Q2rInput, templatetexts[3])
    qe = settings["qe"]
    if qe["manager"] == "local"
        bin = qe["bin"]
        manager = LocalManager(qe["n"], true)
    elseif qe["manager"] == "docker"
        n = qe["n"]
        bin = qe["bin"]
        # manager = DockerEnvironment(n, qe["container"], bin)
    else
    end
    return (
        template = template,
        pressures = settings["pressures"] .* u"GPa",
        dirs = map(settings["pressures"]) do pressure
            abspath(joinpath(
                expanduser(settings["dir"]),
                template[1].control.prefix,
                "p" * string(pressure),
            ))
        end,
        bin = bin,
        manager = manager,
    )
end # function _expand_settings

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

parsecell(str) =
    tryparsefinal(CellParametersCard, str), tryparsefinal(AtomicPositionsCard, str)

end
