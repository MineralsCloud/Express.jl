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

## Loading settings

See
"[Loading settings](https://mineralscloud.github.io/AbInitioSoftwareBase.jl/dev/troubleshooting/#Loading-settings)"
for detailed information.

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
