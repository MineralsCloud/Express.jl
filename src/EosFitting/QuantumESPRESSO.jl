module QuantumESPRESSO

using Crystallography: cellvolume
using Dates: now
using Distributed: LocalManager
using EquationsOfState.Collections
using QuantumESPRESSO.Inputs: inputstring, getoption
using QuantumESPRESSO.Inputs.PWscf: CellParametersCard, PWInput, optconvert
using QuantumESPRESSO.Outputs.PWscf:
    Preamble, parse_electrons_energies, parsefinal, isjobdone
using Setfield: @set!
using Unitful: NoUnits, @u_str, ustrip
using UnitfulAtomic: bohr, Ry

using ...Express:
    Step,
    SelfConsistentField,
    VariableCellOptimization,
    Prepare,
    Analyse,
    _uparse,
    calculationtype

import ...Express
import ..EosFitting

EosFitting.getpotentials(template::PWInput) =
    [x.pseudopot for x in template.atomic_species.data]

EosFitting.getpotentialdir(template::PWInput) = expanduser(template.control.pseudo_dir)

function EosFitting._set_press_vol(template::PWInput, pressure, volume)
    @set! template.cell.press = ustrip(u"kbar", pressure)
    factor = cbrt(volume / (cellvolume(template) * bohr^3)) |> NoUnits  # This is dimensionless and `cbrt` works with units.
    if template.cell_parameters === nothing || getoption(template.cell_parameters) == "alat"
        @set! template.system.celldm[1] *= factor
    else
        @set! template.system.celldm = zeros(6)
        @set! template.cell_parameters =
            optconvert("bohr", CellParametersCard(template.cell_parameters.data * factor))
    end
    return template
end # function EosFitting.set_press_vol

function EosFitting._check_software_settings(settings)
    map(("manager", "bin", "n")) do key
        @assert haskey(settings, key) "key `$key` not found!"
    end
    @assert isinteger(settings["n"]) && settings["n"] >= 1
    if settings["manager"] == "docker"
        @assert haskey(settings, "container")
    elseif settings["manager"] == "ssh"
    elseif settings["manager"] == "local"  # Do nothing
    else
        error("unknown manager `$(settings["manager"])`!")
    end
end # function _check_software_settings

const EosMap = (
    m = Murnaghan,
    bm2 = BirchMurnaghan2nd,
    bm3 = BirchMurnaghan3rd,
    bm4 = BirchMurnaghan4th,
    v = Vinet,
)

function Express.Settings(settings)
    template = parse(PWInput, read(expanduser(settings["template"]), String))
    qe = settings["qe"]
    if qe["manager"] == "local"
        n = qe["n"]
        bin = qe["bin"]
        manager = LocalManager(n, true)
    elseif qe["manager"] == "docker"
        n = qe["n"]
        bin = qe["bin"]
        # manager = DockerEnvironment(n, qe["container"], bin)
    else
    end
    return (
        template = template,
        pressures = settings["pressures"] .* u"GPa",
        trial_eos = EosMap[Symbol(settings["trial_eos"]["type"])](settings["trial_eos"]["parameters"] .*
                                                                  _uparse.(settings["trial_eos"]["units"])...),
        dirs = map(settings["pressures"]) do pressure
            abspath(joinpath(
                expanduser(settings["dir"]),
                template.control.prefix,
                "p" * string(pressure),
            ))
        end,
        manager = manager,
    )
end # function Settings

function EosFitting.preset(step, template, args...)
    @set! template.control.verbosity = "high"
    @set! template.control.wf_collect = true
    @set! template.control.tstress = true
    @set! template.control.tprnfor = true
    @set! template.control.disk_io = "high"
    @set! template.control.calculation =
        calculationtype(step) <: SelfConsistentField ? "scf" : "vc-relax"
    @set! template.control.outdir = join(
        [
            template.control.prefix,
            template.control.calculation,
            string(now()),
            string(rand(UInt)),
        ],
        "_",
    )
    return template
end

function (::EosFitting.Step{SelfConsistentField,Analyse{:output}})(s::AbstractString)
    return parse(Preamble, s).omega => parse_electrons_energies(s, :converged).ε[end]  # volume, energy
end
function (::EosFitting.Step{VariableCellOptimization,Analyse{:output}})(s::AbstractString)
    if !isjobdone(s)
        @warn "Job is not finished!"
    end
    return cellvolume(parsefinal(CellParametersCard{Float64}, s)) =>
        parse_electrons_energies(s, :converged).ε[end]  # volume, energy
end

Express.inputstring(object::PWInput) = inputstring(object)

end # module QuantumESPRESSO
