module QuantumESPRESSO

using AbInitioSoftwareBase.Inputs: set_verbosity
using Crystallography: Cell, eachatom, cellvolume
using Distributed: LocalManager
using EquationsOfState.Collections
using QuantumESPRESSO.Inputs: inputstring, optionof
using QuantumESPRESSO.Inputs.PWscf:
    CellParametersCard, AtomicPositionsCard, PWInput, optconvert
using QuantumESPRESSO.Outputs.PWscf:
    Preamble, parse_electrons_energies, parsefinal, isjobdone, tryparsefinal
using QuantumESPRESSO.CLI: PWCmd
using Setfield: @set!
using Unitful
using UnitfulAtomic

using ...Express: SelfConsistentField, VariableCellOptimization
import ..EosFitting:
    preset_template,
    _check_software_settings,
    _expand_settings,
    _readoutput

export safe_exit

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
                                                                  uparse.(
            settings["trial_eos"]["units"];
            unit_context = [Unitful, UnitfulAtomic],
        )...),
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

function preset_template(calc, template)
    template = set_verbosity(template, "high")
    @set! template.control.calculation = calc isa SelfConsistentField ? "scf" : "vc-relax"
    @set! template.control.outdir = mktempdir()
    return template
end

function _readoutput(::SelfConsistentField, s::AbstractString)
    preamble = tryparse(Preamble, s)
    e = try
        parse_electrons_energies(s, :converged)
    catch
        nothing
    end
    if preamble !== nothing && e !== nothing
        return preamble.omega * u"bohr^3" => e.ε[end] * u"Ry"  # volume, energy
    else
        return
    end
end # function _readoutput
function _readoutput(::VariableCellOptimization, s::AbstractString)
    if !isjobdone(s)
        @warn "Job is not finished!"
    end
    x = tryparsefinal(CellParametersCard, s)
    if x !== nothing
        return cellvolume(parsefinal(CellParametersCard, s)) * u"bohr^3" =>
            parse_electrons_energies(s, :converged).ε[end] * u"Ry"  # volume, energy
    else
        return
    end
end

safe_exit(template::PWInput, dir) = touch(joinpath(dir, template.control.prefix * ".EXIT"))

end # module QuantumESPRESSO
