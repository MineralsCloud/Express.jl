"""
# module StructureOptimization



# Examples

```jldoctest
julia>
```
"""
module StructureOptimization

using LinearAlgebra
using Statistics

using IntervalArithmetic
using IntervalRootFinding
using EquationsOfState
using EquationsOfState
using QuantumESPRESSOBase
using QuantumESPRESSOBase.QuantumESPRESSOInput.PW
using QuantumESPRESSOParsers.InputParsers
using Setfield
using Slurm.Commands
using Slurm.Scriptify
using Slurm.Shells

export update_alat,
    generate_input,
    generate_script

const ALAT_LENS = @lens _.system.celldm[1]

function update_alat(pw::PWInput, alats::AbstractVecOrMat)
    results = similar(alats, typeof(pw))
    for (i, alat) in enumerate(alats)
        results[i] = set(pw, ALAT_LENS, alat)
    end
    return results
end
function update_alat(pw::PWInput, eos::EquationOfState, pressures::AbstractVecOrMat)
    volumes = [find_volume(eos, p, 0..1000, Newton).interval.lo for p in pressures]
    alats = (volumes ./ det(pw.cellparameters.data)).^(1 / 3)
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
    for (out, input) in zip(output, update_alat(pw, eos, pressures))
        write(out, input, debug)
    end
end

function generate_script(shell::Shell, sbatch::Sbatch, modules, pressures::AbstractVecOrMat)
    content = "#!$(string(shell.path))\n"
    content *= scriptify(sbatch) * '\n'
    content *= scriptify(modules) * '\n'
    content *= """
    for i in $(join(pressures, ' ')); do
        (
            cd vc_\$i || return
            mpirun -np $(div(24, length(pressures))) pw.x -in vc_\$i.in > vc_\$i.out &
            sleep 3
        )
    done
    wait
    """
    open(sbatch.script, "r+") do io
        write(io, content)
    end
end # function generate_script

function total(pw::PWInput, trial_eos::EquationOfState, pressures::AbstractVecOrMat)
    eos = fit_energy(trial_eos, volumes, parse_total_energy(outfiles))
    generate_input(pw, eos, pressures, new_input)
end # function total

end
