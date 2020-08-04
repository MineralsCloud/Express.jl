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
import ..Phonon: preset, prep_input, preprocess, _expand_settings

function prep_input(::DfptMethod, template::PhInput, from::PWInput)
    template = preset_template(template)
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
function preset_template(template::PWInput)
    @set! template.control.calculation = "scf"
    return set_verbosity(template, "high")
end # function preset
function preset_template(template::PhInput)
    @set! template.inputph.verbosity = "high"
    return template
end # function preset

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
