# Installation

## Install Julia

First, you should install Julia. We recommend downloading it from
[its official website](https://julialang.org/downloads/). Please follow the detailed instructions on its website if you have to
[build Julia from source](https://github.com/JuliaLang/julia/blob/master/doc/build/build.md).
Some computing centers provide preinstalled Julia. Please contact your administrator for
more information in that case.

If you have [Homebrew](https://brew.sh) installed, open
`Terminal.app` and type

```shell
$ brew install --cask julia
```

if you want to install it as a prebuilt binary app, or

```shell
$ brew install julia
```

if you want to install it as a [formula](https://docs.brew.sh/Formula-Cookbook).

If you want to install multiple Julia versions in the same operating system,
a suggested way is to use a version manager such as
[`asdf`](https://asdf-vm.com/guide/introduction.html).
First, [install `asdf`](https://asdf-vm.com/guide/getting-started.html#_3-install-asdf).
Then, run

```shell
$ asdf install julia 1.6.5   # or other versions of Julia
$ asdf global julia 1.6.5
```

to install Julia and
[set `v1.6.5` as a global version](https://asdf-vm.com/guide/getting-started.html#_6-set-a-version).

You can also try another cross-platform installer for the Julia programming language [`juliaup`](https://github.com/JuliaLang/juliaup).

### Compatibility

Versions higher than `v1.3`,
especially `v1.6`, are strongly recommended. This package may not work on `v1.0` and below.
Since the Julia team has set `v1.6` as the long-term support (LTS) release, support for versions below `v1.6` will gradually drop.

Julia and Julia packages support multiple operating systems and CPU architectures; check [this table](https://julialang.org/downloads/#supported_platforms) to see if it can be installed on your machine. For Mac computers with M-series processors, this package and its dependencies may not work. Please install the Intel-compatible version of Julia (for macOS x86).

## Install `Express`

Now I am using [macOS](https://en.wikipedia.org/wiki/MacOS) as a standard
platform to explain the following steps:

1. Open `Terminal.app`, and type `julia` to start an interactive session (known as
   [REPL](https://docs.julialang.org/en/v1/stdlib/REPL/)).

2. Run the following commands and wait for them to finish:

   ```@repl
   using Pkg; Pkg.update()
   Pkg.add("Express")
   ```

3. Run

   ```@repl
   using Express
   ```

   and have fun!

4. While using, please keep this Julia session alive. Restarting might recompile
   the package and cost some time.

If you want to install the latest in development (maybe buggy) version of `Express`, type

```@repl
using Pkg; Pkg.update()
pkg"add Express#master"
```

in the second step instead.

## Uninstall and reinstall `Express`

1. To uninstall, in a Julia session, run

   ```@repl
   using Pkg; Pkg.rm("Express"); Pkg.gc()
   ```

2. Press `ctrl+d` to quit the current session. Start a new Julia session and
   reinstall `Express`.
