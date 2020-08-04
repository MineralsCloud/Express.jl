module QuantumESPRESSO

using AbInitioSoftwareBase.Inputs: inputstring, writeinput, set_verbosity
using Distributed: LocalManager
using QuantumESPRESSO.Inputs.PWscf:
    AtomicPositionsCard, CellParametersCard, PWInput, optconvert
using QuantumESPRESSO.Inputs.PHonon: PhInput, Q2rInput, MatdynInput, DynmatInput, relayinfo
using QuantumESPRESSO.Outputs.PWscf: tryparsefinal
using Setfield: @set!
using Unitful: @u_str
using UnitfulAtomic

using ...Express: DfptMethod, SelfConsistentField, ForceConstant
import ..Phonon: preset_template, _expand_settings, parsecell

# This is a helper function and should not be exported.
function preset_template(::DfptMethod, template::PhInput, pw::PWInput)
    @set! template.inputph.verbosity = "high"
    return relayinfo(pw, template)
end
function preset_template(::SelfConsistentField, template::PWInput)
    @set! template.control.calculation = "scf"
    @set! template.control.outdir = mktempdir()
    return set_verbosity(template, "high")
end

# function (::Step{ForceConstant,Action{:prepare_input}})(
#     q2r_input,
#     phonon_input,
#     template::Q2rInput,
# )
#     object = parse(PhInput, read(phonon_input, String))
#     write(q2r_input, inputstring(relay(object, template)))
#     return
# end

function _expand_settings(settings)
    templatetexts = [read(expanduser(f), String) for f in settings["template"]]
    template = parse(PWInput, templatetexts[1]),
    parse(PhInput, templatetexts[2]),
    parse(Q2rInput, templatetexts[3])
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

parsecell(str) =
    tryparsefinal(CellParametersCard, str), tryparsefinal(AtomicPositionsCard, str)

end
