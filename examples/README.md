# How to run this example

Install the latest version of Julia (as new as you can).
Install [Docker Desktop](https://docs.docker.com/docker-for-mac/install/), start it.

Start your Julia terminal, run

```julia
julia> using Pkg

julia> pkg"add https://github.com/MineralsCloud/Pseudopotentials.jl.git"

julia> pkg"add https://github.com/singularitti/LibSymspg.jl.git"

julia> pkg"add https://github.com/MineralsCloud/Crystallography.jl.git"

julia> pkg"add https://github.com/singularitti/DockerPy.jl.git"

julia> pkg"add https://github.com/singularitti/PyFortran90Namelists.jl.git"

julia> pkg"add https://github.com/MineralsCloud/QuantumESPRESSOBase.jl.git"

julia> pkg"add https://github.com/MineralsCloud/QuantumESPRESSOParsers.jl.git"

julia> pkg"add https://github.com/MineralsCloud/QuantumESPRESSO.jl.git"

julia> pkg"dev https://github.com/MineralsCloud/Express.jl.git"
```

After the last command, this repo will be cloned to `DEPOT_PATH`. On a *nix system, it is
usually `~/.julia/dev/`. Go to this `~/.julia/dev/Express`, start your VSCode or Atom at
**exactly** there. Open a Julia REPL in your VSCode/Atom, run the content in
`examples/GaN.jl` line by line. When it is calling Quantum ESPRESSO to do the calculations,
it looks like the REPL stucks, but it is not. Just wait! Do not stop the REPL! Each scf
calculation will take about 2 minutes on 2 processors and each vc-relax calculation will
take about 6-9 minutes. So it might need 2-3 hours to run the whole workflow, depending on
how good your computer is.

## Troubleshooting

If you run into a `ConnectionError(..., PermissionError(..., 'Permission denied'))` reported
by some `.py` files on a Linux system like Ubuntu. A possible solution is to follow
[this tutorial](https://docs.docker.com/engine/install/linux-postinstall/#manage-docker-as-a-non-root-user).
(A restart of your system may be required.) After that you do not need root privileges to
run `docker` commands, and try the same commands again.
