"""
# module StructureOptimization



# Examples

```jldoctest
julia>
```
"""
module StructureOptimization

using LinearAlgebra

using EquationsOfState
using QuantumESPRESSOBase
using QuantumESPRESSOBase.QuantumESPRESSOInput.PW
using QuantumESPRESSOParsers.InputParsers
using Setfield

export update_alat,
    generate_input

const ALAT_LENS = @lens _.system.celldm[1]

function update_alat(pw::PWInput, alats::AbstractVecOrMat)
    results = similar(alats, typeof(pw))
    for (i, alat) in enumerate(alats)
        results[i] = set(pw, ALAT_LENS, alat)
    end
    return results
end
function update_alat(pw::PWInput, eos::EquationOfState, pressures::AbstractVecOrMat)
    volumes = find_volume(eos, pressures)
    alats = (volumes ./ det(pw.cellparameters)).^(1 / 3)
    update_alat(pw, alats)
end

function generate_input(
    pw::PWInput,
    eos::EquationOfState,
    pressures::AbstractVecOrMat,
    output::AbstractVecOrMat,
    debug::Bool = true
)
    @assert size(output) == size(pressures)
    for input in update_alat(pw, eos, pressures)
        write(output, input, debug)
    end
end

function total(pw::PWInput, trial_eos::EquationOfState, pressures::AbstractVecOrMat)
    eos = fit_energy(trial_eos, volumes, parse_total_energy(outfiles))
    generate_input(pw, eos, pressures, new_input)
end # function total

end
