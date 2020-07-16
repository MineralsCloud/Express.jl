# How to run this example

Install the latest version of Julia (as new as you can). Install this package as instructed
in ["Installation"](@ref) section.

After the last command, this repo will be cloned to `DEPOT_PATH`. On a *nix system, it is
usually `~/.julia/dev/`. Go to this `~/.julia/dev/Express`, start your VS Code or Atom at
**exactly** there. Open a Julia REPL in your VS Code/Atom, run

```julia
julia> using Express.EosFitting

julia> preprocess(SelfConsistentField(), "<PATH-TO-EXPRESS>/examples/GaN/eos.yaml")

julia> jobs = process(SelfConsistentField(), "<PATH-TO-EXPRESS>/examples/GaN/eos.yaml")
```

Then it will start running scf calculations on different pressures. You can type `jobs`
in terminal at any time to inquire the status of the jobs. When they are running,
a `🚧` emoji will be shown. If succeeded, a `✅` will be shown. If failed, a `❌` will
be shown.

When it is calling Quantum ESPRESSO to do the calculations,
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
