module QuantumESPRESSO

using Crystallography: Cell, eachatom, cellvolume
using Dates: now
using Distributed: LocalManager
using EquationsOfState.Collections
using OptionalArgChecks: @argcheck
using QuantumESPRESSO.Inputs: inputstring, getoption
using QuantumESPRESSO.Inputs.PWscf:
    CellParametersCard,
    AtomicPositionsCard,
    PWInput,
    optconvert,
    set_verbosity,
    set_structure
using QuantumESPRESSO.Outputs.PWscf:
    Preamble, parse_electrons_energies, parsefinal, isjobdone, tryparsefinal
using Setfield: @set!
using Unitful: NoUnits, @u_str, ustrip
using UnitfulAtomic

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

export safe_exit

EosFitting.getpotentials(template::PWInput) =
    [x.pseudopot for x in template.atomic_species.data]

EosFitting.getpotentialdir(template::PWInput) = expanduser(template.control.pseudo_dir)

function EosFitting._set_press_vol(template::PWInput, pressure, volume)::PWInput
    @set! template.cell.press = ustrip(u"kbar", pressure)
    factor = cbrt(volume / (cellvolume(template) * u"bohr^3")) |> NoUnits  # This is dimensionless and `cbrt` works with units.
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
        @argcheck haskey(settings, key) "key `$key` not found!"
    end
    @argcheck isinteger(settings["n"]) && settings["n"] >= 1
    if settings["manager"] == "docker"
        @argcheck haskey(settings, "container")
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
        trial_eos = EosMap[Symbol(settings["trial_eos"]["type"])](settings["trial_eos"]["parameters"] .*
                                                                  _uparse.(settings["trial_eos"]["units"])...),
        dirs = map(settings["pressures"]) do pressure
            abspath(joinpath(
                expanduser(settings["dir"]),
                template.control.prefix,
                "p" * string(pressure),
            ))
        end,
        bin = bin,
        manager = manager,
    )
end # function Settings

function EosFitting.preset(step, template, args...)
    template = set_verbosity(template, "high")
    @set! template.control.calculation =
        calculationtype(step) <: SelfConsistentField ? "scf" : "vc-relax"
    @set! template.control.outdir = join(
        [template.control.prefix, template.control.calculation, now(), rand(UInt)],
        "_",
    )
    return template
end

function EosFitting.analyse(step, s::AbstractString)
    if calculationtype(step) <: SelfConsistentField
        return parse(Preamble, s).omega * u"bohr^3" =>
            parse_electrons_energies(s, :converged).ε[end] * u"Ry"  # volume, energy
    else
        if !isjobdone(s)
            @warn "Job is not finished!"
        end
        return cellvolume(parsefinal(CellParametersCard, s)) * u"bohr^3" =>
            parse_electrons_energies(s, :converged).ε[end] * u"Ry"  # volume, energy
    end
end

safe_exit(template::PWInput, dir) = touch(joinpath(dir, template.control.prefix * ".EXIT"))

EosFitting.parsecell(str::AbstractString) =
    tryparsefinal(CellParametersCard, str), tryparsefinal(AtomicPositionsCard, str)

EosFitting.set_structure(template::PWInput, c, a) = set_structure(template, c, a)

end # module QuantumESPRESSO
