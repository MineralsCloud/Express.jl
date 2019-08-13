"""
# module StructureOptimization



# Examples

```jldoctest
julia>
```
"""
module StructureOptimization

using LinearAlgebra

using Kaleido: @batchlens

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
    generate_input!,
    generate_script

function update_alat(pw::PWscfInput, eos::EquationOfState, pressure::Real)
    volume = find_volume(PressureTarget, eos, pressure, 0..1000, Newton).interval.lo
    alat = cbrt(volume / det(pw.cell_parameters.data))

    lenses = @batchlens begin
        _.system.celldm[1]
        _.cell.press
    end
    set(pw, lenses, (alat, pressure))
end

function generate_input!(
    inputs::AbstractVecOrMat,
    template::PWscfInput,
    eos::EquationOfState,
    pressures::AbstractVecOrMat,
    verbose::Bool = true
)
    @assert size(inputs) == size(pressures)
    for (input, objects) in zip(inputs, map(p -> update_alat(template, eos, p), pressures))
        write(input, objects, verbose)
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
