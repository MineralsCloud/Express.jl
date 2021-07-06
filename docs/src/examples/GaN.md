# hcp-GaN example

## Fitting equations of state

### Run interactively

Install the latest version of Julia (as new as you can). Install this package as instructed
in "[Installation](@ref)" section.

After the last command, this repo will be cloned to `DEPOT_PATH`. On a *nix system, it is
usually `~/.julia/dev/`. Go to this `~/.julia/dev/Express`, start your VS Code or Atom at
**exactly** there. Open a Julia REPL in your VS Code/Atom, run

```julia
julia> using Express.EquationOfStateWorkflows

julia> config = "<PATH-TO-EXPRESS>/examples/GaN/eos.yaml";

julia> preprocess(SelfConsistentField(), config)  # Step 1

julia> scfjobs = process(SelfConsistentField(), config)  # Step 2
```

Then it will start running scf calculations on different pressures. You can type `scfjobs`
in terminal at any time to inquire the status of the jobs. When they are running,
a `🚧` emoji will be shown. If succeeded, a `✅` will be shown. If failed, a `❌` will
be shown.

When all jobs are finished (you need at least 6 pressures finished), run

```julia
julia> postprocess(SelfConsistentField(), config)  # Step 3
EquationsOfState.Collections.BirchMurnaghan3rd{Unitful.Quantity{Float64,D,U} where U where D}
 v0 = 317.7479585715598 a₀³
 b0 = 172.85031797933803 GPa
 b′0 = 4.379649372725796
 e0 = -612.4315064611411 Ry
```

An equation of state, with units, will be shown.

The next 3 steps are basically the same, just with `VariableCellOptimization` as calculation
type.

```julia
julia> preprocess(VariableCellOptimization(), config)  # Step 1

julia> vcjobs = process(VariableCellOptimization(), config)  # Step 2

julia> postprocess(VariableCellOptimization(), config)
EquationsOfState.Collections.BirchMurnaghan3rd{Unitful.Quantity{Float64,D,U} where U where D}
 v0 = 317.7517433287492 a₀³
 b0 = 172.74897940782353 GPa
 b′0 = 4.388034458575274
 e0 = -612.4315074654779 Ry
```

When it is calling Quantum ESPRESSO to do the calculations,
it looks like the REPL got stuck, but it is not. Just wait! Do not stop the REPL! Each scf
calculation will take about 2 minutes on 2 processors, and each vc-relax calculation will
take about 6-9 minutes. So it might need 2-3 hours to run the whole workflow, depending on
how good your computer is.

### Run using a script

When running on a high performance computer, it is always recommended to first start
a debugging/developing environment, do a small test to make sure that inputs can be
generated and job can be submitted. Then to scale your calculation, write the above commands
in a Julia script and submit it to a schedule manager. That is, write

```julia
using Express.EquationOfStateWorkflows

config = "<PATH-TO-EXPRESS>/examples/GaN/eos.yaml"
preprocess(SelfConsistentField(), config)  # Step 1
scfjobs = process(SelfConsistentField(), config)  # Step 2
print(postprocess(SelfConsistentField(), config))  # Step 3
preprocess(VariableCellOptimization(), config)  # Step 4
vcjobs = process(VariableCellOptimization(), config)  # Step 5
print(postprocess(VariableCellOptimization(), config))  # Step 6
```

to a `.jl` file, say `GaN.jl`.
If you are using [Slurm](https://slurm.schedmd.com/documentation.html),
then write a `job.sh` with

```sh
#!/bin/sh
#SBATCH -A <your-account>
#SBATCH -N 1
#SBATCH --tasks-per-node=24
#SBATCH -J job
#SBATCH --time=02:00:00

julia <PATH-TO-GaN.jl>
```

Then submit it using `sbatch job.sh`.
