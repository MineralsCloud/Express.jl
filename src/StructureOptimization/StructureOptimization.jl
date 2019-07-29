"""
# module StructureOptimization



# Examples

```jldoctest
julia>
```
"""
module StructureOptimization

using EquationsOfState
using QuantumESPRESSO.QuantumESPRESSOInput.PW
using Setfield

export update_alat

const ALAT_LENS = @lens _.system.celldm[1]

function update_alat(pw::PWInput, alats::AbstractVecOrMat)
    results = similar(alats, typeof(pw))
    for (i, alat) in enumerate(alats)
        results[i] = set(pw, ALAT_LENS, alat)
    end
    return results
end
function update_alat(pw::PWInput, eos::EquationOfState, pressures::AbstractVecOrMat)
    volumes = eval_volume(eos, pressures)
    alats = (volumes ./ det(pw.cellparameters)).^(1 / 3)
    update_alat(pw, alats)
end

end
