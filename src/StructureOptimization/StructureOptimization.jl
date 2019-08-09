"""
# module StructureOptimization



# Examples

```jldoctest
julia>
```
"""
module StructureOptimization

using LinearAlgebra

using IntervalArithmetic
using IntervalRootFinding
using EquationsOfState
using EquationsOfState.Collections
using EquationsOfState.NumericallyFindVolume
using QuantumESPRESSOBase
using QuantumESPRESSOBase.Inputs.PWscf
using QuantumESPRESSOParsers.InputParsers
using Setfield
using SlurmWorkloadFileGenerator.Commands
using SlurmWorkloadFileGenerator.Scriptify
using SlurmWorkloadFileGenerator.Shells

export update_alat,
    generate_input,
    generate_script

function update_alat(pw::PWscfInput, eos::EquationOfState, pressure::Real)
    volume = find_volume(PressureTarget, eos, pressure, 0..1000, Newton).interval.lo
    alat = cbrt(volume / det(pw.cell_parameters.data))

    alat_lens = @lens _.system.celldm[1]
    press_lens = @lens _.cell.press
    set(set(pw, alat_lens, alat), press_lens, pressure)
end

function generate_input(
    pw::PWscfInput,
    eos::EquationOfState,
    pressures::AbstractVecOrMat,
    output::AbstractVecOrMat,
    debug::Bool = true
)
    @assert size(output) == size(pressures)
    for (out, objects) in zip(output, map(v -> update_alat(pw, eos, v), pressures))
        write(out, objects, debug)
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

function total(pw::PWscfInput, trial_eos::EquationOfState, pressures::AbstractVecOrMat)
    eos = fit_energy(trial_eos, volumes, parse_total_energy(outfiles))
    generate_input(pw, eos, pressures, new_input)
end # function total

end
