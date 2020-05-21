"""
# module EosFitting



# Examples

```jldoctest
julia>
```
"""
module EosFitting

using Compat: isnothing, only
using ConstructionBase: setproperties, constructorof
using Crystallography
using Crystallography.Arithmetics: cellvolume
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
using JSON
using LinearAlgebra: det
using Parameters: @with_kw
using QuantumESPRESSO.Inputs: InputFile, getoption, qestring
using QuantumESPRESSO.Inputs.PWscf:
    AtomicPositionsCard, CellParametersCard, PWInput, optconvert
using QuantumESPRESSO.Outputs: OutputFile
using QuantumESPRESSO.Outputs.PWscf:
    Preamble, parse_electrons_energies, parsefinal, isjobdone
using QuantumESPRESSO.CLI: PWExec
using QuantumESPRESSO.Setters: VerbositySetter, CellParametersSetter
using Setfield: set, @set!
using Unitful
using UnitfulAtomic
using YAML

using Express:
    ScfCalculation,
    PhononCalculation,
    StructureOptimization,
    CPMD,
    PrepareInput,
    LaunchJob,
    AnalyseOutput

import ..Step
using ..CLI: MpiExec
using ..Jobs: nprocs_task, distribute_process

export Step,
    ScfCalculation,
    StructureOptimization,
    PrepareInput,
    LaunchJob,
    AnalyseOutput,
    InputFile,
    load_settings,
    parse_template,
    set_alat_press

Step(::ScfCalculation, ::PrepareInput) = Step(1)
Step(::ScfCalculation, ::LaunchJob) = Step(2)
Step(::ScfCalculation, ::AnalyseOutput) = Step(3)
Step(::StructureOptimization, ::PrepareInput) = Step(4)
Step(::StructureOptimization, ::LaunchJob) = Step(5)
Step(::StructureOptimization, ::AnalyseOutput) = Step(6)

function _check_settings(settings)
    map(("template", "np", "pressures", "trial_eos", "qe", "workdir")) do key
        @assert haskey(settings, key) "`$key` is reuqired but not found in settings!"
    end
    if length(settings["qe"]) > 1
        error("multiple Quantum ESPRESSO methods are given! It must be 1!")
    end
    @assert only(keys(settings["qe"])) ∈ ("local", "docker", "ssh")
    @assert isdir(settings["workdir"])
    @assert isfile(settings["template"])
    @assert isinteger(settings["np"]) && settings["np"] >= 1
    if length(settings["pressures"]) <= 6
        @info "pressures less than 6 may give unreliable results, consider more if possible!"
    end
    map(("type", "parameters")) do key
        @assert haskey(settings["trial_eos"], key) "`$key` is reuqired for eos `$type`!"
    end
end # function _check_settings

const EosMap = (
    m = Murnaghan,
    bm2 = BirchMurnaghan2nd,
    bm3 = BirchMurnaghan3rd,
    bm4 = BirchMurnaghan4th,
)

function _expand_settings(settings)
    template = parse_template(InputFile(abspath(expanduser(settings["template"]))))
    return (
        template = template,
        pressures = settings["pressures"] * u"GPa",
        trial_eos = getindex(EosMap, Symbol(settings["trial_eos"]["type"]))(settings["trial_eos"]["parameters"]...),
        inputs = map(settings["pressures"]) do pressure
            abspath(joinpath(
                expanduser(settings["workdir"]),
                "p" * string(pressure),
                template.control.calculation,
                template.control.prefix * ".in",
            ))
        end,
        np = settings["np"],
        qe = settings["qe"]
    )
end # function _expand_settings

function load_settings(path::AbstractString)
    settings = _loadfrom(path)
    _check_settings(settings)  # Errors will be thrown if exist
    return _expand_settings(settings)
end # function load_settings

parse_template(str::AbstractString) = parse(PWInput, str)
parse_template(file::InputFile) = parse(PWInput, read(file))

function set_alat_press(
    template::PWInput,
    eos::EquationOfState{<:Unitful.AbstractQuantity},
    pressure::Unitful.AbstractQuantity,
)
    volume = findvolume(eos(Pressure()), pressure, (eps(float(eos.v0)), 1.3 * eos.v0))  # In case `eos.v0` has a `Int` as `T`. See https://github.com/PainterQubits/Unitful.jl/issues/274.
    factor = cbrt(volume / (cellvolume(template) * u"bohr^3")) |> NoUnits  # This is dimensionless and `cbrt` works with units.
    if isnothing(template.cell_parameters)
        @set! template.system.celldm[1] *= factor
    else
        @set! template.system.celldm = zeros(6)
        @set! template.cell_parameters =
            optconvert("bohr", template.cell_parameters * factor)
    end
    @set! template.cell.press = ustrip(u"kbar", pressure)
    return template
end # function update_alat_press

# This is a helper function and should not be exported.
_preset(step::Union{Step{1},Step{4}}, template::PWInput) = setproperties(
    template,
    control = setproperties(
        template.control,
        calculation = step isa Step{1} ? "scf" : "vc-relax",
        verbosity = "high",
        tstress = true,
        tprnfor = true,
    ),
)

function (step::Union{Step{1},Step{4}})(
    inputs,
    template::PWInput,
    trial_eos::EquationOfState,
    pressures,
)
    template = _preset(step, template)
    map(inputs, pressures) do input, pressure  # `map` will check size mismatch
        mkpath(joinpath(splitpath(input)[1:end-1]...))
        touch(input)
        object = set_alat_press(template, trial_eos, pressure)  # Create a new `object` from `template`, with its `alat` and `pressure` changed
        write(InputFile(input), object)  # Write the `object` to a Quantum ESPRESSO input file
        object
    end
end
function (::Step{1})(path::AbstractString)
    template, pressures, trial_eos, inputs = load_settings(path)
    return Step(1)(inputs, template, trial_eos, pressures)
end # function preprocess

function _dockercmd(exec::PWExec, n, input)
    "sh -c 'mpiexec --mca btl_vader_single_copy_mechanism none -np $n " *
    string(exec())[2:end-1] *
    " -inp $input'"
end # function _wrapcmd

function (::Union{Step{2},Step{5}})(
    inputs,
    outputs,
    np,
    exec::PWExec,
    ids = workers();
    isdocker::Bool = true,
    container = nothing,
)
    # `map` guarantees they are of the same size, no need to check.
    n = nprocs_task(np, length(inputs))
    cmds = map(inputs, outputs) do input, output  # A vector of `Cmd`s
        if isdocker
            _dockercmd(exec, n, input)
        else
            MpiExec(n)(exec(stdin = input, stdout = output))
        end
    end
    return distribute_process(
        cmds,
        ids;
        isdocker = isdocker,
        container = container,
        inputs = inputs,
    )
end
function (step::Union{Step{2},Step{5}})(path::AbstractString)
    settings = load_settings(path)
    outputs = map(Base.Fix2(replace, ".in" => ".out"), settings.inputs)
    return step(settings.inputs, outputs, settings.np, PWExec(which = settings.qe["bin"]))
end

function (step::Union{Step{3},Step{6}})(
    outputs,
    trial_eos::EquationOfState{<:Unitful.AbstractQuantity},
)
    xy = map(outputs) do output
        s = read(OutputFile(output))
        if !isjobdone(s)
            @warn "Job is not finished!"
        end
        _results(step, s), parse_electrons_energies(s, :converged).ε[end]  # volume, energy
    end
    return lsqfit(trial_eos(Energy()), first.(xy) .* u"bohr^3", last.(xy) .* u"Ry")
end # function postprocess

_results(::Step{3}, s::AbstractString) = parse(Preamble, s).omega
_results(::Step{6}, s::AbstractString) =
    cellvolume(parsefinal(CellParametersCard{Float64}, s))

function _saveto(filepath::AbstractString, data)
    ext = _extension(filepath)
    if ext ∈ (".yaml", ".yml")
        YAML.write_file(expanduser(filepath), data)
    elseif ext == ".json"
        open(expanduser(filepath), "w") do io
            JSON.print(io, data)
        end
    else
        error("unknown file extension `$ext`!")
    end
end # function _saveto

function _loadfrom(filepath::AbstractString)
    ext = _extension(filepath)
    if ext ∈ (".yaml", ".yml")
        return open(expanduser(filepath), "r") do io
            YAML.load(io)
        end
    elseif ext == ".json"
        return JSON.parsefile(expanduser(filepath))
    else
        error("unknown file extension `$ext`!")
    end
end # function _loadfrom

_extension(filepath::AbstractString) = filepath |> splitext |> last |> lowercase

_uparse(str::AbstractString) = uparse(str; unit_context = [Unitful, UnitfulAtomic])

end
