module QuantumESPRESSO

using Crystallography: Cell, eachatom, cellvolume
using Dates: now
using Distributed: LocalManager
using EquationsOfState.Collections
using QuantumESPRESSO.Inputs: inputstring, getoption
using QuantumESPRESSO.Inputs.PWscf:
    CellParametersCard,
    AtomicPositionsCard,
    PWInput,
    optconvert,
    set_verbosity,
    set_structure,
    set_press_vol
using QuantumESPRESSO.Outputs.PWscf:
    Preamble, parse_electrons_energies, parsefinal, isjobdone, tryparsefinal
using QuantumESPRESSO.CLI: PWCmd
using Setfield: @set!
using Unitful: @u_str
using UnitfulAtomic

using ...Express: SelfConsistentField, VariableCellOptimization, _uparse
using ..EosFitting: set_press_vol
import ..EosFitting:
    getpotentials,
    getpotentialdir,
    _set_press_vol,
    _set_structure,
    _check_software_settings,
    _expand_settings,
    _prep_input,
    analyse,
    parsecell,
    set_structure

export safe_exit

getpotentials(template::PWInput) = [x.pseudopot for x in template.atomic_species.data]

getpotentialdir(template::PWInput) = expanduser(template.control.pseudo_dir)

_set_press_vol(template::PWInput, pressure, volume) =
    set_press_vol(template, pressure, volume)

_set_structure(template::PWInput, cell_parameters, atomic_positions) =
    set_structure(template, cell_parameters, atomic_positions)

function _check_software_settings(settings)
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

function _expand_settings(settings)
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
        bin = PWCmd(; bin = bin),
        manager = manager,
    )
end # function _expand_settings

function _prep_input(calculation, template)
    template = set_verbosity(template, "high")
    @set! template.control.calculation =
        calculation isa SelfConsistentField ? "scf" : "vc-relax"
    @set! template.control.outdir = join(
        [template.control.prefix, template.control.calculation, now(), rand(UInt)],
        "_",
    )
    return template
end

analyse(::SelfConsistentField, s::AbstractString) =
    parse(Preamble, s).omega * u"bohr^3" =>
        parse_electrons_energies(s, :converged).ε[end] * u"Ry"  # volume, energy
function analyse(::VariableCellOptimization, s::AbstractString)
    if !isjobdone(s)
        @warn "Job is not finished!"
    end
    return cellvolume(parsefinal(CellParametersCard, s)) * u"bohr^3" =>
        parse_electrons_energies(s, :converged).ε[end] * u"Ry"  # volume, energy
end

safe_exit(template::PWInput, dir) = touch(joinpath(dir, template.control.prefix * ".EXIT"))

parsecell(str) =
    tryparsefinal(CellParametersCard, str), tryparsefinal(AtomicPositionsCard, str)

end # module QuantumESPRESSO
