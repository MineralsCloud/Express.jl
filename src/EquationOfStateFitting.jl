"""
# module EquationOfStateFitting



# Examples

```jldoctest
julia>
```
"""
module EquationOfStateFitting

using LinearAlgebra: det

using Distributed: workers
using Compat: isnothing
using EquationsOfState
using EquationsOfState.Collections: EquationOfState
using EquationsOfState.NonlinearFitting: lsqfit
using EquationsOfState.Find: findvolume
using Kaleido: @batchlens
using QuantumESPRESSO: to_qe, cell_volume
using QuantumESPRESSO.Cards: optionof
using QuantumESPRESSO.Cards.PWscf: AtomicPositionsCard, CellParametersCard
using QuantumESPRESSO.Inputs: autofill_cell_parameters
using QuantumESPRESSO.Inputs.PWscf: PWInput
using QuantumESPRESSO.Outputs.PWscf:
    Preamble, parse_electrons_energies, parsefinal, isjobdone
using QuantumESPRESSOBase.CLI: PWCmd
using Setfield: set, @lens
using Unitful
using UnitfulAtomic

import ..Step
using ..Jobs: MpiCmd, nprocs_per_subjob, distribute_process

export update_alat_press, preprocess, postprocess, submit

function update_alat_press(
    template::PWInput,
    eos::EquationOfState{<:Unitful.AbstractQuantity},
    pressure::Unitful.AbstractQuantity,
)
    if isnothing(template.cell_parameters)
        template = autofill_cell_parameters(template)
    end
    # In case `eos.v0` has a `Int` as `T`. See https://github.com/PainterQubits/Unitful.jl/issues/274.
    v0 = float(eos.v0)
    volume = findvolume(PressureForm(), eos, pressure, (eps(v0), 1.3v0))
    option = optionof(template.cell_parameters)
    # Please do not change this logic, see https://github.com/MineralsCloud/Express.jl/pull/38.
    determinant = if option ∈ ("alat", "bohr")
        # `alat` uses relative values WRT `celldm`, which uses "bohr" as unit.
        # So `"alat"` is equivalent to `"bohr"`.
        det(template.cell_parameters.data) * u"bohr^3"
    elseif option == "angstrom"
        det(template.cell_parameters.data) * u"angstrom^3" |> u"bohr^3"
    else  # This should never happen actually.
        error("unknown option $option in a `CellParametersCard`!")
    end
    # `cbrt` works with units.
    alat = cbrt(volume / determinant) |> NoUnits  # This is dimensionless.
    lenses = @batchlens(begin
        _.system.celldm  # Get the `template`'s `system.celldm` value
        _.cell.press     # Get the `template`'s `cell.press` value
        _.cell_parameters.option
    end)
    # Set the `template`'s `system.celldm` and `cell.press` values
    return set(template, lenses, ([alat], ustrip(u"kbar", pressure), "alat"))
end # function update_alat_press

# This is a helper function and should not be exported.
function _boilerplate(step::Step{N}, template::PWInput) where {N}
    lenses = @batchlens(begin
        _.control.calculation  # Get the `template`'s `control.calculation` value
        _.control.verbosity    # Get the `template`'s `control.verbosity` value
        _.control.tstress      # Get the `template`'s `control.tstress` value
        _.control.tprnfor      # Get the `template`'s `control.tprnfor` value
    end)
    # Set the `template`'s values with...
    return set(template, lenses, (N == 1 ? "scf" : "vc-relax", "high", true, true))
end # function _boilerplate

function preprocess(
    step::Step,
    inputs::AbstractArray{<:AbstractString},
    template::PWInput,
    trial_eos::EquationOfState,
    pressures::AbstractArray,
    verbose::Bool = false,
)
    @assert(
        size(inputs) == size(pressures),
        "The `inputs` and `pressures` must be of the same size!"
    )  # `zip` does not guarantee they are of the same size, must check explicitly.
    template = _boilerplate(step, template)
    for (input, pressure) in zip(inputs, pressures)
        # Get a new `object` from the `template`, with its `alat` and `pressure` changed
        object = update_alat_press(template, trial_eos, pressure)
        # `write` will create a file if it doesn't exist.
        write(input, to_qe(object, verbose = verbose))  # Write the `object` to a Quantum ESPRESSO input file
    end
    return
end # function preprocess

function submit(
    inputs::AbstractVector{<:AbstractString},
    outputs::AbstractVector{<:AbstractString},
    np::Int,
    cmdtemplate::MpiCmd,
    ids::AbstractVector{<:Integer} = workers(),
)
    each = nprocs_per_subjob(np, length(inputs))
    if isnothing(cmdtemplate)
        cmdtemplate = MpiCmd(np = each, subcmd = PWCmd(inp = inputs[1]))
    end
    cmds = fill(cmdtemplate, size(inputs))
    for (i, (input, output)) in enumerate(zip(inputs, outputs))
        cmds[i] = set(cmds[i], @lens(_.np), each)
        cmds[i] = set(cmds[i], @lens(_.subcmd), PWCmd(inp = input))
        cmds[i] = set(cmds[i], @lens(_.stdout), output)
    end
    return distribute_process(cmds, ids)
end # function submit

function postprocess(
    ::Step{N},
    outputs::AbstractVector{<:AbstractString},
    trial_eos::EquationOfState{<:Unitful.AbstractQuantity},
) where {N}
    energies, volumes = zeros(length(outputs)), zeros(length(outputs))
    for (i, output) in enumerate(outputs)
        open(output, "r") do io
            s = read(io, String)
            isjobdone(s) || @warn("Job is not finished!")
            energies[i] = parse_electrons_energies(s, :converged).ε[end]
            volumes[i] = if N == 1
                parse(Preamble, s).omega
            elseif N == 2
                cell_volume(parsefinal(CellParametersCard, s))
            else
                error("The step $N must be `1` or `2`!")
            end
        end
    end
    return lsqfit(EnergyForm(), trial_eos, volumes .* u"bohr^3", energies .* u"Ry")
end # function postprocess

end
