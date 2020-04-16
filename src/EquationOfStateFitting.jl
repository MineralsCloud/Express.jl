"""
# module EquationOfStateFitting



# Examples

```jldoctest
julia>
```
"""
module EquationOfStateFitting

using Compat: isnothing
using Crystallography.Arithmetics: cellvolume
using Distributed: workers
using EquationsOfState.Collections: EquationOfState, Pressure, Energy
using EquationsOfState.NonlinearFitting: lsqfit
using EquationsOfState.Find: findvolume
using Kaleido: @batchlens
using LinearAlgebra: det
using QuantumESPRESSO.Inputs: InputFile, getoption, qestring
using QuantumESPRESSO.Inputs.PWscf: AtomicPositionsCard, CellParametersCard, PWInput
using QuantumESPRESSO.Outputs: OutputFile
using QuantumESPRESSO.Outputs.PWscf: Preamble, parse_electrons_energies, parsefinal, isjobdone
using QuantumESPRESSOBase.CLI: PWCmd
using QuantumESPRESSOBase.Setters: VerbositySetter, CellParametersSetter
using Setfield: set, @set!
using Unitful
using UnitfulAtomic

import ..Step
using ..CLI: MpiExec
using ..Jobs: nprocs_task, distribute_process

export update_alat_press, preprocess, postprocess, fire

function update_alat_press(
    template::PWInput,
    eos::EquationOfState{<:Unitful.AbstractQuantity},
    pressure::Unitful.AbstractQuantity,
)
    if isnothing(template.cell_parameters)
        template = set(template, CellParametersSetter())
    end
    volume = findvolume(eos(Pressure()), pressure, (eps(float(eos.v0)), 1.3 * eos.v0))  # In case `eos.v0` has a `Int` as `T`. See https://github.com/PainterQubits/Unitful.jl/issues/274.
    alat = cbrt(volume / cellvolume(template) * u"bohr^3") |> NoUnits  # This is dimensionless and `cbrt` works with units.
    lenses = @batchlens(begin
        _.system.celldm  # Get the `template`'s `system.celldm` value
        _.cell.press     # Get the `template`'s `cell.press` value
        _.cell_parameters.option
    end)
    return set(template, lenses, ([alat], ustrip(u"kbar", pressure), "alat"))  # Set the `template`'s `system.celldm` and `cell.press` values
end # function update_alat_press

# This is a helper function and should not be exported.
function _preset(step::Step{N}, template::PWInput) where {N}
    lenses = @batchlens(begin
        _.control.calculation  # Get the `template`'s `control.calculation` value
        _.control.verbosity    # Get the `template`'s `control.verbosity` value
        _.control.tstress      # Get the `template`'s `control.tstress` value
        _.control.tprnfor      # Get the `template`'s `control.tprnfor` value
    end)
    # Set the `template`'s values with...
    return set(template, lenses, (N == 1 ? "scf" : "vc-relax", "high", true, true))
end # function _preset

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
    end  # `zip` does not guarantee they are of the same size, must check explicitly.
    template = _preset(step, template)
    objects = similar(inputs, PWInput)  # Create an array of `undef` of `PWInput` type
    for (i, (input, pressure)) in enumerate(zip(inputs, pressures))
        # Create a new `object` from the `template`, with its `alat` and `pressure` changed
        object = update_alat_press(template, trial_eos, pressure)
        # `write` will create a file if it doesn't exist.
        objects[i] = object
        write(InputFile(input), object)  # Write the `object` to a Quantum ESPRESSO input file
    end
    return objects
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
    end  # `zip` does not guarantee they are of the same size, must check explicitly.
    n = nprocs_task(np, length(inputs))
    cmds = similar(inputs, Base.AbstractCmd)
    for (i, (input, output)) in enumerate(zip(inputs, outputs))
        lenses = @batchlens(begin
            _.n
            _.cmd.inp
        end)
        cmds[i] = pipeline(convert(Cmd, set(template, lenses, (n, input))), stdout = output)
    end
    return distribute_process(cmds, ids)
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
