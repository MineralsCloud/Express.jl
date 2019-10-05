"""
# module EquationOfStateFitting



# Examples

```jldoctest
julia>
```
"""
module EquationOfStateFitting

using LinearAlgebra: det

using Compat: isnothing
using EquationsOfState
using EquationsOfState.Collections: EquationOfState
using EquationsOfState.NonlinearFitting: lsqfit
using EquationsOfState.Find: findvolume
using Kaleido: @batchlens
using MLStyle: @match
using QuantumESPRESSOBase: to_qe, option
using QuantumESPRESSOBase.Inputs.PWscf: PWscfInput, autofill_cell_parameters
using QuantumESPRESSOParsers.OutputParsers.PWscf: parse_total_energy,
                                                  parse_cell_parameters,
                                                  parse_head,
                                                  isjobdone
using Setfield: set
using Unitful: AbstractQuantity, ustrip, NoUnits, @u_str
using UnitfulAstro
using UnitfulAtomic

import ..Step
using ..SelfConsistentField: write_metadata

export update_alat_press, prepare, finish

function update_alat_press(
    template::PWscfInput,
    eos::EquationOfState,
    pressure::AbstractQuantity,
)
    # In case `eos.v0` has a `Int` as `T`. See https://github.com/PainterQubits/Unitful.jl/issues/274.
    isnothing(template.cell_parameters) && (template = autofill_cell_parameters(template))
    v0 = float(eos.v0)
    volume = findvolume(PressureForm(), eos, pressure, (eps(v0), 1.3v0))
    # If the `CellParametersCard` contains a matrix of plain numbers (no unit).
    determinant = @match option(template.cell_parameters) begin
        # `alat` uses relative values WRT `celldm`, which uses "bohr" as unit.
        # So `"alat"` is equivalent to `"bohr"`.
        "alat" || "bohr" => det(cell_parameters.data) * u"bohr^3"
        "angstrom" => det(template.cell_parameters.data) * u"angstrom^3" |> u"bohr^3"
    end
    # `cbrt` works with units.
    alat = cbrt(volume / determinant) |> NoUnits  # This is dimensionless.
    lenses = @batchlens(begin
        _.system.celldm ∘ _[$1]  # Get the `template`'s `system.celldm[1]` value
        _.cell.press             # Get the `template`'s `cell.press` value
        _.cell_parameters.option
    end)
    # Set the `template`'s `system.celldm[1]` and `cell.press` values with `alat` and `pressure`
    return set(template, lenses, (alat, ustrip(u"kbar", pressure), "alat"))
end # function update_alat_press

# This is a helper function and should not be exported.
function _preset(step::Step{N}, template::PWscfInput) where {N}
    lenses = @batchlens(begin
        _.control.calculation  # Get the `template`'s `control.calculation` value
        _.control.verbosity    # Get the `template`'s `control.verbosity` value
        _.control.tstress      # Get the `template`'s `control.tstress` value
        _.control.tprnfor      # Get the `template`'s `control.tprnfor` value
    end)
    # Set the `template`'s values with...
    return set(template, lenses, (N == 1 ? "scf" : "vc-relax", "high", true, true))
end # function _preset

function prepare(
    step::Step,
    inputs::AbstractVector{<:AbstractString},
    template::PWscfInput,
    trial_eos::EquationOfState,
    pressures::AbstractVector,
    metadatafiles::AbstractVector{<:AbstractString} = map(x -> splitext(x)[1] * ".json", inputs),
    verbose::Bool = false
)
    # Check parameters
    @assert(
        length(inputs) == length(pressures) == length(metadatafiles),
        "The inputs, pressures and the metadata files must be the same size!"
    )
    template = _preset(step, template)
    # Write input and metadata
    for (input, pressure, metadatafile) in zip(inputs, pressures, metadatafiles)
        # Get a new `object` from the `template`, with its `alat` and `pressure` changed
        object = update_alat_press(template, trial_eos, pressure)
        write(input, to_qe(object, verbose = verbose))  # Write the `object` to a Quantum ESPRESSO input file
        write_metadata(metadatafile, input, template)
    end
    return
end # function prepare

function finish(
    ::Step{N},
    outputs::AbstractVector{<:AbstractString},
    trial_eos::EquationOfState,
) where {N}
    energies, volumes = zeros(length(outputs)), zeros(length(outputs))
    for (i, output) in enumerate(outputs)
        open(output, "r") do io
            s = read(io, String)
            isjobdone(s) || @warn("Job is not finished!")
            energies[i] = parse_total_energy(s)[end]
            volumes[i] = @match N begin
                1 => parse_head(s)["unit-cell volume"]
                2 => parse_cell_parameters(s)[end] |> det
                _ => error("The step $N must be `1` or `2`!")
            end
        end
    end
    return lsqfit(EnergyForm(), trial_eos, volumes, energies)
end # function finish

# This is what the web service does
# function workflow(args)
#     prepare(Step(1), inputs, template, trial_eos, pressures, metadatafiles)
#     trial_eos = finish(outputs, trial_eos)
#     prepare(Step(2), inputs, template, trial_eos, pressures, metadatafiles)
#     return finish(outputs, trial_eos)
# end # function workflow

end
