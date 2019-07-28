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

end
