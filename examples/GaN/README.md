# How to run this example

Install the latest version of Julia (as new as you can).
Install [Docker Desktop](https://docs.docker.com/docker-for-mac/install/), start it.

Start your Julia terminal, run

```julia
using Pkg
pkg"add git@github.com:MineralsCloud/AbInitioSoftwareBase.jl.git"
pkg"add git@github.com:MineralsCloud/Pseudopotentials.jl.git"
pkg"add git@github.com:MineralsCloud/Crystallography.jl.git"
pkg"add git@github.com:singularitti/DockerPy.jl.git"
pkg"add git@github.com:singularitti/PyFortran90Namelists.jl.git"
pkg"add git@github.com:MineralsCloud/QuantumESPRESSOBase.jl.git"
pkg"add git@github.com:MineralsCloud/QuantumESPRESSOParsers.jl.git"
pkg"add git@github.com:MineralsCloud/QuantumESPRESSO.jl.git"
pkg"dev git@github.com:MineralsCloud/Express.jl.git"
```

After the last command, this repo will be cloned to `DEPOT_PATH`. On a *nix system, it is
usually `~/.julia/dev/`. Go to this `~/.julia/dev/Express`, start your VS Code or Atom at
**exactly** there. Open a Julia REPL in your VS Code/Atom, run the content in
`examples/GaN.jl` line by line. When it is calling Quantum ESPRESSO to do the calculations,
it looks like the REPL got stuck, but it is not. Just wait! Do not stop the REPL! Each scf
calculation will take about 2 minutes on 2 processors, and each vc-relax calculation will
take about 6-9 minutes. So it might need 2-3 hours to run the whole workflow, depending on
how good your computer is.

## Troubleshooting

If you run into a `ConnectionError(..., PermissionError(..., 'Permission denied'))` reported
by some `.py` files on a Linux system like Ubuntu, a possible solution is to follow
[this tutorial](https://docs.docker.com/engine/install/linux-postinstall/#manage-docker-as-a-non-root-user).
(A restart of your system may be required.) After that, you do not need root privileges to
run `docker` commands, and try the same commands again.
