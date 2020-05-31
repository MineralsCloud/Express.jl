"""
# module EosFitting



# Examples

```jldoctest
julia>
```
"""
module EosFitting

using Crystallography: cellvolume
using EquationsOfState.Collections: Pressure, Energy, EquationOfState
using EquationsOfState.NonlinearFitting: lsqfit
using EquationsOfState.Find: findvolume
using Unitful: NoUnits
using UnitfulAtomic: bohr, Ry
using QuantumESPRESSO.CLI: pwcmd
using QuantumESPRESSO.Inputs: InputFile, inputstring

using ..Express:
    Step,
    SelfConsistentField,
    VariableCellRelaxation,
    PrepareInput,
    LaunchJob,
    AnalyseOutput,
    load_settings,
    _uparse
using ..Jobs: nprocs_task, distribute_process
using ..Workspaces: DockerWorkspace, LocalWorkspace

export Step,
    SelfConsistentField,
    VariableCellRelaxation,
    PrepareInput,
    LaunchJob,
    AnalyseOutput,
    InputFile,
    load_settings,
    parse_template,
    set_press_vol

function set_press_vol(template, pressure, eos; minscale = eps(), maxscale = 1.3)
    @assert minscale > zero(minscale)  # No negative volume
    volume = findvolume(eos(Pressure()), pressure, (minscale, maxscale) .* eos.v0)
    return set_press_vol(template, pressure, volume)
end # function set_press_vol

function (step::Step{T,PrepareInput})(
    inputs,
    template,
    pressures,
    trial_eos::EquationOfState;
    kwargs...,
) where {T}
    template = _set_verbosity(T, template)
    map(inputs, pressures) do input, pressure  # `map` will check size mismatch
        mkpath(dirname(input))
        object = set_press_vol(template, pressure, trial_eos; kwargs...)  # Create a new `object` from `template`, with its `alat` and `pressure` changed
        open(input, "w") do io
            write(io, inputstring(object))
        end
    end
end
function (step::Step{T,PrepareInput})(path::AbstractString) where {T}
    settings = load_settings(path)
    return step(settings.inputs, settings.template, settings.pressures, settings.trial_eos)
end # function preprocess

function (::Step{T,LaunchJob})(outputs, inputs, workspace, cmdtemplate) where {T}
    # `map` guarantees they are of the same size, no need to check.
    n = nprocs_task(workspace.n, length(inputs))
    cmds = map(inputs, outputs) do input, output  # A vector of `Cmd`s
        _generate_cmds(n, cmdtemplate, input, workspace)
    end
    return distribute_process(cmds, workspace)
end
function (step::Step{T,LaunchJob})(path::AbstractString) where {T}
    settings = load_settings(path)
    outputs = map(Base.Fix2(replace, ".in" => ".out"), settings.inputs)
    return step(settings.inputs, outputs, settings.nprocs, pwcmd(bin = settings.qe["bin"]))
end

function (step::Step{T,AnalyseOutput})(outputs, trial_eos) where {T}
    xy = map(outputs) do output
        s = read(output, String)
        parseenergies(step, s)
    end
    return lsqfit(trial_eos(Energy()), first.(xy) .* bohr^3, last.(xy) .* Ry)
end # function postprocess

# function (::T)(
#     outputs,
#     inputs,
#     template,
#     pressures,
#     trial_eos,
#     workspace,
#     cmd,
# ) where {T<:Union{SelfConsistentField,VariableCellRelaxation}}
#     Step{typeof(T),PrepareInput}(inputs, template, pressures, trial_eos)
#     Step{typeof(T),LaunchJob}(outputs, inputs, workspace, cmd)
#     Step{typeof(T),AnalyseOutput}(outputs, trial_eos)
# end

_generate_cmds(n, cmd, input, ::DockerWorkspace) =
    "sh -c 'mpiexec --mca btl_vader_single_copy_mechanism none -np $n " *
    string(cmd)[2:end-1] *
    " -inp $input'"

function parse_template end

function _results end

function tostring end

function parseenergies end

function _set_verbosity end

module QuantumESPRESSO

using Compat: isnothing
using Crystallography: cellvolume
using ConstructionBase: setproperties
using EquationsOfState.Collections
using QuantumESPRESSO.Inputs: InputFile, inputstring, getoption
using QuantumESPRESSO.Inputs.PWscf: CellParametersCard, PWInput, optconvert
using QuantumESPRESSO.Outputs: OutputFile
using QuantumESPRESSO.Outputs.PWscf:
    Preamble, parse_electrons_energies, parsefinal, isjobdone
using QuantumESPRESSO.CLI: pwcmd
using Setfield: @set!
using Unitful: NoUnits, @u_str, ustrip
using UnitfulAtomic: bohr

using ...Express:
    Step,
    SelfConsistentField,
    VariableCellRelaxation,
    PrepareInput,
    LaunchJob,
    AnalyseOutput,
    load_settings,
    _uparse
using ...Jobs: nprocs_task, distribute_process
using ...Workspaces: DockerWorkspace

import ...Express
import ..EosFitting

function EosFitting.set_press_vol(template::PWInput, pressure, volume)
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

function _check_qe_settings(settings)
    map(("scheme", "bin")) do key
        @assert haskey(settings, key)
    end
    if settings["scheme"] == "docker"
        @assert haskey(settings, "container")
    elseif settings["scheme"] == "ssh"
    elseif settings["scheme"] == "local"  # Do nothing
    else
        error("unknown scheme `$(settings["scheme"])`!")
    end
end # function _check_qe_settings

function Express._check_settings(settings)
    map(("template", "nprocs", "pressures", "trial_eos", "qe", "dir")) do key
        @assert haskey(settings, key)
    end
    _check_qe_settings(settings["qe"])
    @assert isdir(settings["dir"])
    @assert isfile(settings["template"])
    @assert isinteger(settings["nprocs"]) && settings["nprocs"] >= 1
    if length(settings["pressures"]) <= 6
        @info "pressures less than 6 may give unreliable results, consider more if possible!"
    end
    map(("type", "parameters", "units")) do key
        @assert haskey(settings["trial_eos"], key)
    end
end # function _check_settings

const EosMap = (
    m = Murnaghan,
    bm2 = BirchMurnaghan2nd,
    bm3 = BirchMurnaghan3rd,
    bm4 = BirchMurnaghan4th,
)

function Express.Settings(settings)
    template = parse_template(InputFile(abspath(expanduser(settings["template"]))))
    return (
        template = template,
        pressures = settings["pressures"] .* u"GPa",
        trial_eos =
            EosMap[Symbol(settings["trial_eos"]["type"])](settings["trial_eos"]["parameters"] .*
                                                          _uparse.(settings["trial_eos"]["units"])...),
        inputs = map(settings["pressures"]) do pressure
            abspath(joinpath(
                expanduser(settings["dir"]),
                "p" * string(pressure),
                template.control.calculation,
                template.control.prefix * ".in",
            ))
        end,
        nprocs = settings["nprocs"],
        qe = settings["qe"],
    )
end # function Settings

EosFitting.parse_template(str::AbstractString) = parse(PWInput, str)
EosFitting.parse_template(file::InputFile) = parse(PWInput, read(file))

# This is a helper function and should not be exported.
function EosFitting._set_verbosity(T, template::PWInput)
    @set! template.control.verbosity = "high"
    @set! template.control.wf_collect = true
    @set! template.control.tstress = true
    @set! template.control.tprnfor = true
    @set! template.control.disk_io = "high"
    @set! template.control.calculation = T <: SelfConsistentField ? "scf" : "vc-relax"
    return template
end # macro _set_verbosity

EosFitting._results(::Step{SelfConsistentField,AnalyseOutput}, s::AbstractString) =
    parse(Preamble, s).omega
EosFitting._results(::Step{VariableCellRelaxation,AnalyseOutput}, s::AbstractString) =
    cellvolume(parsefinal(CellParametersCard{Float64}, s))

function EosFitting.parseenergies(step, s)
    if !isjobdone(s)
        @warn "Job is not finished!"
    end
    return _results(step, s), parse_electrons_energies(s, :converged).ε[end]  # volume, energy
end # function EosFitting.parseenergies

end # module QuantumESPRESSO

end
