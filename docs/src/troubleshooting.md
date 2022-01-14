# Troubleshooting

```@contents
Pages = ["troubleshooting.md"]
```

This page collects some possible errors you may encounter and trick how to fix them.

*If you have additional tips, please submit a PR with suggestions.*

## Installation problems

### Cannot find the Julia executable

Make sure you have Julia installed in your environment. Please download the latest
[stable version](https://julialang.org/downloads/) for your platform.
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

### Have trouble installing [`PyCall.jl`](https://github.com/JuliaPy/PyCall.jl)

If you are seeing error message like this:

```
ERROR: LoadError: InitError: PyError (PyImport_ImportModule

The Python package qha could not be imported by pyimport. Usually this means
that you did not install qha in the Python version being used by PyCall.
```

One solution is to re-configure `PyCall` to use a different Python
version on your system: set `ENV["PYTHON"]` to the path of the python
executable you want to use, run `Pkg.build("PyCall")`, and re-launch Julia.
Please see [this part](https://github.com/JuliaPy/PyCall.jl#specifying-the-python-version)
for more detailed instructions.

Another solution is to configure `PyCall` to use a Julia-specific Python
distribution via the [`Conda.jl` package](https://github.com/JuliaPy/Conda.jl)
(which installs a private Anaconda
Python distribution), which has the advantage that packages can be installed
and kept up-to-date via Julia.  As explained in the `PyCall` documentation, in Julia,
set `ENV["PYTHON"]=""`, run `Pkg.build("PyCall")`, and re-launch Julia.

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
Otherwise check other `YAML` syntax you may have broken.

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

Then you probably forget loading a plugin package for `Express.jl`. For example, you should
run

```julia-repl
julia> using QuantumESPRESSOExpress
```

to fix this error.
