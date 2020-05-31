"""
# module EosFitting



# Examples

```jldoctest
julia>
```
"""
module EosFitting

using Compat: isnothing
using ConstructionBase: setproperties
using Crystallography: cellvolume
using Distributed: workers
using EquationsOfState.Collections:
    EquationOfState,
    Pressure,
    Energy,
    BirchMurnaghan2nd,
    BirchMurnaghan3rd,
    BirchMurnaghan4th,
    Murnaghan
using EquationsOfState.NonlinearFitting: lsqfit
using EquationsOfState.Find: findvolume
using QuantumESPRESSO.Inputs: InputFile
using QuantumESPRESSOBase.Inputs.PWscf: CellParametersCard, PWInput, optconvert
using QuantumESPRESSO.Outputs: OutputFile
using QuantumESPRESSO.Outputs.PWscf:
    Preamble, parse_electrons_energies, parsefinal, isjobdone
using QuantumESPRESSOBase.CLI: pwcmd
using Setfield: set, @set!
using Unitful
using UnitfulAtomic

using ..Express:
    Step,
    SelfConsistentField,
    VariableCellRelaxation,
    PrepareInput,
    LaunchJob,
    AnalyseOutput,
    load,
    save,
    _uparse
using ..CLI: mpicmd
using ..Jobs: nprocs_task, distribute_process
using ..Workspaces: DockerWorkspace

export Step,
    SelfConsistentField,
    VariableCellRelaxation,
    PrepareInput,
    LaunchJob,
    AnalyseOutput,
    InputFile,
    load_settings,
    parse_template,
    set_vol_press

function set_vol_press(template::PWInput, eos, pressure)
    volume = findvolume(eos(Pressure()), pressure, (eps(float(eos.v0)), 1.3 * eos.v0))  # In case `eos.v0` has a `Int` as `T`. See https://github.com/PainterQubits/Unitful.jl/issues/274.
    factor = cbrt(volume / (cellvolume(template) * u"bohr^3")) |> NoUnits  # This is dimensionless and `cbrt` works with units.
    if isnothing(template.cell_parameters)
        @set! template.system.celldm[1] *= factor
    else
        if template.cell_parameters.option == "alat"
            @set! template.system.celldm[1] *= factor
        else
            @set! template.system.celldm = zeros(6)
            @set! template.cell_parameters = optconvert(
                "bohr",
                CellParametersCard(template.cell_parameters.data * factor),
            )
        end
    end
    @set! template.cell.press = ustrip(u"kbar", pressure)
    return template
end # function set_vol_press

function (step::Step{T,PrepareInput})(inputs, template, trial_eos, pressures) where {T}
    template = _preset(step, template)
    map(inputs, pressures) do input, pressure  # `map` will check size mismatch
        mkpath(joinpath(splitpath(input)[1:end-1]...))
        touch(input)
        object = set_vol_press(template, trial_eos, pressure)  # Create a new `object` from `template`, with its `alat` and `pressure` changed
        write(InputFile(input), object)  # Write the `object` to a Quantum ESPRESSO input file
    end
end
function (step::Step{T,PrepareInput})(path::AbstractString) where {T}
    template, pressures, trial_eos, inputs = load_settings(path)
    return step(inputs, template, trial_eos, pressures)
end # function preprocess

function (::Step{T,LaunchJob})(inputs, outputs, np, cmd, ids = workers()) where {T}
    # `map` guarantees they are of the same size, no need to check.
    n = nprocs_task(np, length(inputs))
    cmds = map(inputs, outputs) do input, output  # A vector of `Cmd`s
    end
    return distribute_process(cmds, ids)
end
function (::Step{T,LaunchJob})(
    inputs,
    outputs,
    workspace::DockerWorkspace,
    cmd,
    ids = workers(),
) where {T}
    # `map` guarantees they are of the same size, no need to check.
    n = nprocs_task(workspace.n, length(inputs))
    cmds = map(inputs, outputs) do input, output  # A vector of `Cmd`s
        _dockercmd(cmd, n, input)
    end
    return distribute_process(workspace, outputs, cmds, ids)
end
function (step::Step{T,LaunchJob})(path::AbstractString) where {T}
    settings = load_settings(path)
    outputs = map(Base.Fix2(replace, ".in" => ".out"), settings.inputs)
    return step(settings.inputs, outputs, settings.nprocs, pwcmd(bin = settings.qe["bin"]))
end

function (step::Step{T,AnalyseOutput})(
    outputs,
    trial_eos::EquationOfState{<:Unitful.AbstractQuantity},
) where {T}
    xy = map(outputs) do output
        s = read(OutputFile(output))
        if !isjobdone(s)
            @warn "Job is not finished!"
        end
        _results(step, s), parse_electrons_energies(s, :converged).ε[end]  # volume, energy
    end
    return lsqfit(trial_eos(Energy()), first.(xy) .* u"bohr^3", last.(xy) .* u"Ry")
end # function postprocess

function (::SelfConsistentField)(
    inputs,
    template,
    trial_eos,
    pressures,
    outputs,
    workspace,
    np,
    cmd,
    ids = workers(),
)
    Step{SelfConsistentField,PrepareInput}(inputs, template, trial_eos, pressures)
    Step{SelfConsistentField,LaunchJob}(outputs, np, cmd, ids)
    Step{SelfConsistentField,AnalyseOutput}(outputs, trial_eos)
end

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

function _check_settings(settings)
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

function Settings(settings)
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

function load_settings(path::AbstractString)
    settings = load(path)
    _check_settings(settings)  # Errors will be thrown if exist
    return Settings(settings)
end # function load_settings

parse_template(str::AbstractString) = parse(PWInput, str)
parse_template(file::InputFile) = parse(PWInput, read(file))

# This is a helper function and should not be exported.
_preset(::Step{T,PrepareInput}, template::PWInput) where {T} = setproperties(
    template,
    control = setproperties(
        template.control,
        calculation = T isa SelfConsistentField ? "scf" : "vc-relax",
        verbosity = "high",
        tstress = true,
        tprnfor = true,
    ),
)

function _dockercmd(exec, n, input)
    "sh -c 'mpiexec --mca btl_vader_single_copy_mechanism none -np $n " *
    string(exec)[2:end-1] *
    " -inp $input'"
end # function _wrapcmd

_results(::Step{SelfConsistentField,AnalyseOutput}, s::AbstractString) =
    parse(Preamble, s).omega
_results(::Step{VariableCellRelaxation,AnalyseOutput}, s::AbstractString) =
    cellvolume(parsefinal(CellParametersCard{Float64}, s))

end
