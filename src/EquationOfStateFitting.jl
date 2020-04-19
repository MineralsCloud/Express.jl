"""
# module EquationOfStateFitting



# Examples

```jldoctest
julia>
```
"""
module EquationOfStateFitting

using Compat: isnothing
using ConstructionBase: setproperties
using Crystallography.Arithmetics: cellvolume
using Distributed: workers
using EquationsOfState.Collections: EquationOfState, Pressure, Energy
using EquationsOfState.NonlinearFitting: lsqfit
using EquationsOfState.Find: findvolume
using JSON
using LinearAlgebra: det
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

export init_settings, load_settings, parse_template, set_alat_press, preprocess, postprocess, fire

function init_settings(path::AbstractString)
    settings = Dict{String,Any}(
        "template" => "",
        "pressures" => zeros(1),
        "trial_eos" => nothing,
        "path" => "",
        "prefix" => "",
    )
    path = expanduser(path)
    ext = path |> splitext |> last |> lowercase
    if ext ∈ (".yaml", ".yml")
        YAML.write_file(path, settings)
    elseif ext == ".json"
        open(path, "w") do io
            JSON.print(io, settings)
        end
    else
        error("unknown file extension `$ext`!")
    end
end # function init_settings

function load_settings(path::AbstractString)
    path = expanduser(path)
    ext = path |> splitext |> last |> lowercase
    if ext ∈ (".yaml", ".yml")
        settings = open(path, "r") do io
            YAML.load(io)
        end
    elseif ext == ".json"
        settings = JSON.parsefile(path)
    else
        error("unknown file extension `$ext`!")
    end
    _check_settings(settings)
    return settings
end # function load_settings

# This is a helper function and should not be exported!
function _check_settings(settings::AbstractDict)
    if isempty(settings["template"])
        @warn "the path of the `template` file is not set!"
    end
    if eltype(settings["pressures"]) ∉ (Int, Float64)
        @warn ""
    end
    if isnothing(settings["trial_eos"])
        @warn "the trial eos is not set!"
    end
    for key in ("path", "prefix")
        if isempty(settings[key])
            @info "key `$key` is not set, will use the default value!"
        end
    end
end # function _check_settings

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
    inputs::AbstractArray{<:AbstractString},
    template::PWInput,
    trial_eos::EquationOfState,
    pressures::AbstractArray,
    verbose::Bool = false,
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

end
