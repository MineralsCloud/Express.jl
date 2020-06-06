module QuantumESPRESSO

using Compat: isnothing
using Crystallography: cellvolume
using ConstructionBase: setproperties
using EquationsOfState.Collections
using QuantumESPRESSO.Inputs: InputFile, inputstring, getoption
using QuantumESPRESSO.Inputs.PWscf: CellParametersCard, PWInput, optconvert
using QuantumESPRESSO.Outputs.PWscf:
    Preamble, parse_electrons_energies, parsefinal, isjobdone
using Setfield: @set!
using Unitful: NoUnits, @u_str, ustrip
using UnitfulAtomic: bohr

using ...Express: Step, SelfConsistentField, VariableCellRelaxation, AnalyseOutput, _uparse
using ...Environments: DockerEnvironment, LocalEnvironment
using ..EosFitting: parse_template

import ...Express
import ..EosFitting

function EosFitting._set_press_vol(template::PWInput, pressure, volume)
    @set! template.cell.press = ustrip(u"kbar", pressure)
    factor = cbrt(volume / (cellvolume(template) * bohr^3)) |> NoUnits  # This is dimensionless and `cbrt` works with units.
    if isnothing(template.cell_parameters) || getoption(template.cell_parameters) == "alat"
        @set! template.system.celldm[1] *= factor
    else
        @set! template.system.celldm = zeros(6)
        @set! template.cell_parameters =
            optconvert("bohr", CellParametersCard(template.cell_parameters.data * factor))
    end
    return template
end # function EosFitting.set_press_vol

function EosFitting._check_software_settings(settings)
    map(("environment", "bin", "n")) do key
        @assert haskey(settings, key) "key `$key` not found!"
    end
    @assert isinteger(settings["n"]) && settings["n"] >= 1
    if settings["environment"] == "docker"
        @assert haskey(settings, "container")
    elseif settings["environment"] == "ssh"
    elseif settings["environment"] == "local"  # Do nothing
    else
        error("unknown environment `$(settings["environment"])`!")
    end
end # function _check_software_settings

const EosMap = (
    m = Murnaghan,
    bm2 = BirchMurnaghan2nd,
    bm3 = BirchMurnaghan3rd,
    bm4 = BirchMurnaghan4th,
)

function Express.Settings(settings)
    template = parse_template(InputFile(abspath(expanduser(settings["template"]))))
    qe = settings["qe"]
    if qe["environment"] == "local"
        n = qe["n"]
        bin = qe["bin"]
        environment = LocalEnvironment(n, bin, ENV)
    elseif qe["environment"] == "docker"
        n = qe["n"]
        bin = qe["bin"]
        environment = DockerEnvironment(n, qe["container"], bin)
    else
    end
    return (
        template = template,
        pressures = settings["pressures"] .* u"GPa",
        trial_eos = EosMap[Symbol(settings["trial_eos"]["type"])](settings["trial_eos"]["parameters"] .*
                                                                  _uparse.(settings["trial_eos"]["units"])...),
        inputs = map(settings["pressures"]) do pressure
            abspath(joinpath(
                expanduser(settings["dir"]),
                "p" * string(pressure),
                template.control.calculation,
                template.control.prefix * ".in",
            ))
        end,
        environment = environment,
    )
end # function Settings

EosFitting.parse_template(str::AbstractString) = parse(PWInput, str)
EosFitting.parse_template(file::InputFile) = parse(PWInput, read(file))

# This is a helper function and should not be exported.
function EosFitting._set_boilerplate(T, template::PWInput)
    @set! template.control.verbosity = "high"
    @set! template.control.wf_collect = true
    @set! template.control.tstress = true
    @set! template.control.tprnfor = true
    @set! template.control.disk_io = "high"
    @set! template.control.calculation = T <: SelfConsistentField ? "scf" : "vc-relax"
    return template
end # macro _set_boilerplate

_results(::Step{SelfConsistentField,AnalyseOutput}, s::AbstractString) =
    parse(Preamble, s).omega
_results(::Step{VariableCellRelaxation,AnalyseOutput}, s::AbstractString) =
    cellvolume(parsefinal(CellParametersCard{Float64}, s))

function EosFitting.parseenergies(step, s)
    if !isjobdone(s)
        @warn "Job is not finished!"
    end
    return _results(step, s), parse_electrons_energies(s, :converged).ε[end]  # volume, energy
end # function EosFitting.parseenergies

Express.inputstring(object::PWInput) = inputstring(object)

end # module QuantumESPRESSO
