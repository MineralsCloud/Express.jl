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
using EquationsOfState.Collections: EquationOfState, Pressure, Energy, BirchMurnaghan3rd
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

export Settings,
    Step,
    ScfCalculation,
    StructureOptimization,
    PrepareInput,
    LaunchJob,
    AnalyseOutput,
    InputFile,
    init_settings,
    load_settings,
    parse_template,
    set_alat_press

Step(::ScfCalculation, ::PrepareInput) = Step(1)
Step(::ScfCalculation, ::LaunchJob) = Step(2)
Step(::ScfCalculation, ::AnalyseOutput) = Step(3)
Step(::StructureOptimization, ::PrepareInput) = Step(4)
Step(::StructureOptimization, ::LaunchJob) = Step(5)
Step(::StructureOptimization, ::AnalyseOutput) = Step(6)

@with_kw struct Settings
    template::String = ""
    pressures::AbstractArray{<:Unitful.AbstractQuantity} = zeros(8) .* u"GPa"
    trial_eos::EquationOfState{<:Unitful.AbstractQuantity} =
        BirchMurnaghan3rd(0.0 * u"angstrom^3", 0.0 * u"GPa", 0.0, 0.0 * u"eV")
    dir::String = "."
    np::Int = 1
end

function _todict(settings::Settings)
    eos = settings.trial_eos
    return Dict(
        "template" => expanduser(settings.template),
        "pressures" => Dict(
            "values" => ustrip.(settings.pressures),
            "unit" => only(unique(unit.(settings.pressures))),
        ),
        "trial_eos" => "",
        "dir" => expanduser(settings.dir),
        "np" => settings.np,
    )
end # function _todict

init_settings(path::AbstractString, settings = Settings()) =
    _saveto(path, _todict(settings))

function load_settings(path::AbstractString)
    settings = _loadfrom(path)
    if _isgood(settings)
        return Settings(
            template = expanduser(settings["template"]),
            pressures = settings["pressures"]["values"] .*
                        _uparse(settings["pressures"]["unit"]),
            trial_eos = eval(Meta.parse(settings["trial_eos"])),
            dir = expanduser(settings["dir"]),
            np = settings["np"],
        )
    else
        @warn "some settings are not good! Check your input!"
        return settings
    end
end # function load_settings

# This is a helper function and should not be exported!
function _isgood(settings)
    if isempty(settings["template"])
        @warn "the path of the `template` file is not set!"
        return false
    end
    if eltype(settings["pressures"]["values"]) ∉ (Int, Float64)
        println(eltype(settings["pressures"]))
        @warn "pressure values must be floats or integers!"
        return false
    end
    if isnothing(settings["trial_eos"])
        @warn "the trial eos is not set!"
        return false
    end
    if isempty(settings["dir"])
        @info "key `\"dir\"` is not set, will use the default value!"
    end
    return true
end # function _isgood

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
        @set! template.cell_parameters = optconvert("bohr", template.cell_parameters * factor)
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
    settings = load_settings(path)
    if settings isa Settings
        pressures = settings.pressures
        template = parse_template(InputFile(settings.template))
        inputs = map(pressures) do pressure
            joinpath(
                settings.dir,
                pressure |> ustrip |> round |> string,
                "scf",
                template.control.prefix * ".in.txt",
            )
        end
        return Step(1)(inputs, template, settings.trial_eos, pressures)
    else
        error("an error setting is given!")
    end
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
    if settings isa Settings
        pressures = settings.pressures
        template = parse_template(InputFile(settings.template))
        inputs = map(pressures) do pressure
            joinpath(
                settings.dir,
                pressure |> ustrip |> round |> string,
                "scf",
                template.control.prefix * ".in.txt",
            )
        end
        outputs = map(Base.Fix2(replace, ".in" => ".out"), inputs)
        # return step(inputs, outputs, MpiExec(settings.np, DockerExec()))
    else
        error("an error setting is given!")
    end
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
_results(::Step{6}, s::AbstractString) = cellvolume(parsefinal(CellParametersCard, s))

function _saveto(filepath::AbstractString, data)
    ext = _getext(filepath)
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
    ext = _getext(filepath)
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

_getext(filepath::AbstractString) = filepath |> splitext |> last |> lowercase

_uparse(str::AbstractString) = uparse(str; unit_context = [Unitful, UnitfulAtomic])

end
