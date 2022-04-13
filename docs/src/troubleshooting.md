# Troubleshooting

This page collects some possible errors you may encounter and trick how to fix them.

_If you have additional tips, please submit a PR with suggestions._

## Installation problems

### Cannot find the Julia executable

Make sure you have Julia installed in your environment. Please download the latest
[stable version](https://julialang.org/downloads/#current_stable_release) for your platform.
If you are using macOS, the recommended way is to use [Homebrew](https://brew.sh).
If you do not want to install Homebrew or you are using other platforms that Julia supports,
download the corresponding binaries. And then create a symbolic link `/usr/local/bin/julia`
to the Julia executable. If `/usr/local/bin/` is not in your `$PATH`, export it to your `$PATH`.
Some clusters, like
[`Habanero`](https://confluence.columbia.edu/confluence/display/rcs/Habanero+HPC+Cluster+User+Documentation),
[`Comet`](https://www.sdsc.edu/support/user_guides/comet.html),
or [`Expanse`](https://www.sdsc.edu/services/hpc/expanse/index.html),
already have Julia installed as a module, you may
just `module load julia` to use it. If not, either install by yourself or contact your
administrator.

### Have trouble installing [`qha`](https://github.com/MineralsCloud/qha)

If you are seeing error message like this:

```
ERROR: LoadError: InitError: PyError (PyImport_ImportModule

The Python package qha could not be imported by pyimport. Usually this means
that you did not install qha in the Python version being used by PyCall.

PyCall is currently configured to use the Python version at:

/usr/bin/python3

and you should use whatever mechanism you usually use (apt-get, pip, conda,
etcetera) to install the Python package containing the qha module.
```

It is because on some operating systems, `python` is already installed, and Julia
selects it as the default binary. But that `python` cannot install third-party Python packages.
So [`qha`](https://github.com/MineralsCloud/qha) cannot be automatically installed
by Julia.

One solution is to re-configure [`PyCall`](https://github.com/JuliaPy/PyCall.jl) to use a different Python
version on your system: set `ENV["PYTHON"]` to the path of the python
executable you want to use, run `Pkg.build("PyCall")`, and re-launch Julia.
For example, in Julia REPL, run

```julia-repl
julia> using Pkg

julia> ENV["PYTHON"] = "... path of the python executable ..."
# ENV["PYTHON"] = raw"C:\Python37-x64\python.exe" # example for Windows, "raw" to not have to escape: "C:\\Python37-x64\\python.exe"
# ENV["PYTHON"] = "/usr/bin/python3.7"            # example for *nix

julia> Pkg.build("PyCall");
```

Please see [this part](https://github.com/JuliaPy/PyCall.jl#specifying-the-python-version)
for more detailed instructions.

Another solution is to configure `PyCall` to use a Julia-specific Python
distribution via the [`Conda.jl` package](https://github.com/JuliaPy/Conda.jl)
(which installs a private Anaconda
Python distribution), which has the advantage that packages can be installed
and kept up-to-date via Julia. As explained in the `PyCall` documentation, in Julia,
run

```julia-repl
julia> using Pkg

julia> ENV["PYTHON"] = "" # empty string

julia> Pkg.build("PyCall");
```

Then re-launch Julia.

If `qha` still cannot be installed, go to the Python binary directory you specified
(for the second solution, go to [`$JULIA_DEPOT_PATH/conda/3/bin`](https://docs.julialang.org/en/v1/manual/environment-variables/#JULIA_DEPOT_PATH))
and run the following command

```shell
$ ./pip install qha
```

!!! note

    At least Python 3.6 and above is required to install qha.
    Please read its [manual](https://mineralscloud.github.io/qha/tutorials/installing.html)
    for more information.

## Loading settings

### Error parsing `YAML` files

If you encounter

```
ERROR: expected '<document start>' but found YAML.BlockMappingStartToken at nothing
```

or

```
ERROR: while scanning a simple key at line n, column 0: could not find expected ':' at line n+1, column 0
```

Check whether you have no space between the `YAML` key and its value like
`key:1` or `key:some text`, etc. To correct, change to `key: 1`, `key: some text`, etc.
Otherwise, check other `YAML` syntax you may have broken.

## Loading `Express`

### Why is Julia compiling/loading modules so slow? What can I do?

First, we recommend you download the latest version of Julia. Usually, the newest version
has the best performance.

If you just want Julia to do a simple task and only once, you could start Julia REPL with

```shell
$ julia --compile=min
```

to minimize compilation or

```shell
$ julia --optimize=0
```

to minimize optimizations, or just use both. Or you could make a system image
and run with

```shell
$ julia --sysimage custom-image.so
```

See [Fredrik Ekre's talk](https://youtu.be/IuwxE3m0_QQ?t=313) for details.

## Running workflows

### How to solve the error `come across ProcessFailedException`?

The reasons for this error are manifold.
One common reason is that express did not find the executables of Quantum ESPRESSO or
`mpiexec` (e.g., their parent directory is not in the `PATH` environment variable) or the
template input files. If so, please specify the paths to these executables in the absolute
path format:

```yaml
cli:
  mpi:
    np: 64
    path: /usr/local/bin/gcc-10.2.0/openmpi/bin/mpiexec
  pw:
    path: /usr/local/bin/gcc-9.2.0/qe-6.5/bin/pw.x
  ph:
    path: /usr/local/bin/gcc-9.2.0/qe-6.5/bin/ph.x
  q2r:
    path: /usr/local/bin/gcc-9.2.0/qe-6.5/bin/q2r.x
  matdyn:
    path: /usr/local/bin/gcc-9.2.0/qe-6.5/bin/matdyn.x
template:
  scf: /home/qe/template.in
  dfpt: /home/qe/ph.in
  q2r: /home/qe/q2r.in
  disp: /home/qe/phdos.in
```

Besides, the paths to the directories in the template files could also be absolute to
avoid unexpected errors:

```fortran
&CONTROL
  verbosity = 'high'
  tstress = .true.
  tprnfor = .true.
  outdir = '/home/qe/temp'
  prefix = 'Ge'
  disk_io = 'high'
  pseudo_dir = '/home/qe/pseudo'
/
```


### Why does `buildworkflow` return the same `Workflow` object if I am starting a new workflow?

This happens when running a workflow (say, a phonon workflow, for example) followed by a
workflow (an equation of state workflow in this case). And after the following Julia
commands:

```julia
using Express.PhononWorkflow.Recipes
using QuantumESPRESSOExpress

wf = buildworkflow("phonon.yaml")
```

you find that `wf` still contains the jobs from the previous workflow, and they are probably
all succeeded (because that is when you want to run the next workflow). The cause is likely
in your configuration file (`phonon.yaml`). You probably have the following setting:

```yaml
save:
  status: status.jls
```

where the `status.jls` already exists and is the serialized `Workflow` object saved from
the previous (equation of state) workflow. `Express.jl` will first load this file if
it exists. To fix this, change the value of the `status` key to another path.

## Miscellaneous errors

1. If you see the following error message

   ```julia-repl
   julia> w = buildworkflow("eos.toml");
   ERROR: MethodError: Cannot `convert` an object of type Dict{String, Any} to an object of type AbInitioSoftwareBase.Commands.CommandConfig
   Closest candidates are:
     convert(::Type{T}, ::Intervals.AnchoredInterval{P, T}) where {P, T} at ~/.julia/packages/Intervals/ua9cq/src/anchoredinterval.jl:181
     convert(::Type{T}, ::Intervals.Interval{T}) where T at ~/.julia/packages/Intervals/ua9cq/src/interval.jl:253
     convert(::Type{T}, ::P) where {T, P<:(Polynomials.AbstractPolynomial{T})} at ~/.julia/packages/Polynomials/WvTSC/src/common.jl:434
     ...
   ```

   Then you probably forget loading a plugin package for `Express.jl`. For example, you should run

   ```julia-repl
   julia> using QuantumESPRESSOExpress
   ```

   to fix this error.
