"""
# module EquationOfStateFitting



# Examples

```jldoctest
julia>
```
"""
module EquationOfStateFitting

using Compat: isnothing, only
using ConstructionBase: setproperties, constructorof
using Crystallography.Arithmetics: cellvolume
using Distributed: workers
using EquationsOfState
using EquationsOfState.Collections: EquationOfState, Pressure, Energy, BirchMurnaghan3rd
using EquationsOfState.NonlinearFitting: lsqfit
using EquationsOfState.Find: findvolume
using JSON
using LinearAlgebra: det
using Parameters: @with_kw
using QuantumESPRESSO.Inputs: InputFile, getoption, qestring
using QuantumESPRESSO.Inputs.PWscf: AtomicPositionsCard, CellParametersCard, PWInput
using QuantumESPRESSO.Outputs: OutputFile
using QuantumESPRESSO.Outputs.PWscf:
    Preamble, parse_electrons_energies, parsefinal, isjobdone
using QuantumESPRESSOBase.CLI: PWCmd
using QuantumESPRESSOBase.Setters: VerbositySetter, CellParametersSetter
using Setfield: set
using Unitful
using UnitfulAtomic
using YAML

import ..Step
using ..CLI: MpiExec
using ..Jobs: nprocs_task, distribute_process

export Settings, init_settings, load_settings, parse_template, set_alat_press, preprocess, postprocess, fire

@with_kw struct Settings
    template::String = ""
    pressures::AbstractArray{<:Unitful.AbstractQuantity} = zeros(8) .* u"GPa"
    trial_eos::EquationOfState{<:Unitful.AbstractQuantity} = BirchMurnaghan3rd(0.0u"angstrom^3", 0.0u"GPa", 0.0, 0.0u"eV")
    dir::String = "."
    prefix::String = "scf"
end

function _todict(settings::Settings)
    eos = settings.trial_eos
    return Dict(
        "template" => expanduser(settings.template),
        "pressures" => Dict("values" => ustrip.(settings.pressures), "unit" => only(unique(unit.(settings.pressures)))),
        "trial_eos" => Dict("type" => constructorof(typeof(eos)), "parameters" => ustrip.(fieldvalues(eos)), "units" => unit.(fieldvalues(eos))),
        "dir" => expanduser(settings.dir),
        "prefix" => settings.prefix
    )
end # function _todict

init_settings(path::AbstractString, settings = Settings()) = _saveto(path, _todict(settings))

function load_settings(path::AbstractString)
    settings = _loadfrom(path)
    if _isgood(settings)
        eosdict = settings["trial_eos"]
        return Settings(
            template = expanduser(settings["template"]),
            pressures = settings["pressures"]["values"] .* uparse(settings["pressures"]["unit"]),
            trial_eos = eval(Meta.parse(eosdict["type"]))((eosdict["parameters"] .* uparse.(eosdict["units"]; unit_context = [Unitful, UnitfulAtomic]))...),
            dir = expanduser(settings["dir"]),
            prefix = settings["prefix"]
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
    for key in ("dir", "prefix")
        if isempty(settings[key])
            @info "key `$key` is not set, will use the default value!"
        end
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
    if isnothing(template.cell_parameters)
        template = set(template, CellParametersSetter())
    end
    volume = findvolume(eos(Pressure()), pressure, (eps(float(eos.v0)), 1.3 * eos.v0))  # In case `eos.v0` has a `Int` as `T`. See https://github.com/PainterQubits/Unitful.jl/issues/274.
    alat = cbrt(volume / (cellvolume(template) * u"bohr^3")) |> NoUnits  # This is dimensionless and `cbrt` works with units.
    return setproperties(
        template,
        system = setproperties(template.system, celldm = [alat]),
        cell = setproperties(template.cell, press = ustrip(u"kbar", pressure)),
        cell_parameters = setproperties(template.cell_parameters, option = "alat"),
    )
end # function update_alat_press

# This is a helper function and should not be exported.
_preset(step::Step{N}, template::PWInput) where {N} = setproperties(
    template.control,
    calculation = N == 1 ? "scf" : "vc-relax",
    verbosity = "high",
    tstress = true,
    tprnfor = true,
)

function preprocess(
    step::Step,
    inputs::AbstractArray,
    template::PWInput,
    trial_eos::EquationOfState,
    pressures::AbstractArray,
)
    if size(inputs) != size(pressures)
        throw(DimensionMismatch("`inputs` and `pressures` must be of the same size!"))
    else
        template = _preset(step, template)
        objects = similar(inputs, PWInput)  # Create an array of `undef` of `PWInput` type
        for (i, (input, pressure)) in enumerate(zip(inputs, pressures))
            object = set_alat_press(template, trial_eos, pressure)  # Create a new `object` from the `template`, with its `alat` and `pressure` changed
            objects[i] = object  # `write` will create a file if it doesn't exist.
            write(InputFile(input), object)  # Write the `object` to a Quantum ESPRESSO input file
        end
        return objects
    end  # `zip` does not guarantee they are of the same size, must check explicitly.
end # function preprocess
function preprocess(::Step{1}, path::AbstractString)
    settings = load_settings(path)
    if settings isa Settings
        pressures = settings.pressures
        inputs = map(pressures) do pressure
            joinpath(settings.dir, string(round(pressure)), "scf", settings.prefix * ".in")
        end
        return preprocess(Step(1), inputs, parse_template(InputFile(settings.template)), settings.trial_eos, pressures)
    else
        error("an error setting is given!")
    end
end # function preprocess

function fire(
    inputs::AbstractArray{<:AbstractString},
    outputs::AbstractArray{<:AbstractString},
    np::Int,
    template::MpiExec = MpiExec(n = 1, cmd = PWCmd(inp = "")),
    ids::AbstractArray{<:Integer} = workers(),
)
    if size(inputs) != size(outputs)
        throw(DimensionMismatch("`inputs` and `outputs` must be of the same size!"))
    else
        n = nprocs_task(np, length(inputs))
        cmds = similar(inputs, Base.AbstractCmd)
        for (i, (input, output)) in enumerate(zip(inputs, outputs))
            cmds[i] = pipeline(
                convert(
                    Cmd,
                    setproperties(
                        template,
                        n = n,
                        input = setproperties(template.cmd, inp = input),
                    ),
                ),
                stdout = output,
            )
        end
        return distribute_process(cmds, ids)
    end  # `zip` does not guarantee they are of the same size, must check explicitly.
end # function fire

function postprocess(
    ::Step{N},
    outputs::AbstractVector{<:AbstractString},
    trial_eos::EquationOfState{<:Unitful.AbstractQuantity},
) where {N}
    energies, volumes = zeros(length(outputs)), zeros(length(outputs))
    for (i, output) in enumerate(outputs)
        s = read(OutputFile(output))
        isjobdone(s) || @warn("Job is not finished!")
        energies[i] = parse_electrons_energies(s, :converged).ε[end]
        volumes[i] = if N == 1
            parse(Preamble, s).omega
        elseif N == 2
            cellvolume(parsefinal(CellParametersCard, s))
        else
            error("The step $N must be `1` or `2`!")
        end
    end
    return lsqfit(trial_eos(Energy()), volumes .* u"bohr^3", energies .* u"Ry")
end # function postprocess

fieldvalues(eos::EquationOfState) = [getfield(eos, i) for i in 1:nfields(eos)]

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

end
