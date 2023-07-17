# [Installation Guide](@id installation)

```@contents
Pages = ["installation.md"]
Depth = 2
```

Here are the installation instructions for package
[ExpressBase](https://github.com/MineralsCloud/ExpressBase.jl).
If you have trouble installing it, please refer to our [Troubleshooting](@ref) page
for more information.

## Install Julia

First, you should install [Julia](https://julialang.org/). We recommend downloading it from
[its official website](https://julialang.org/downloads/). Please follow the detailed
instructions on its website if you have to
[build Julia from source](https://docs.julialang.org/en/v1/devdocs/build/build/).
Some computing centers provide preinstalled Julia. Please contact your administrator for
more information in that case.
Here's some additional information on
[how to set up Julia on HPC systems](https://github.com/hlrs-tasc/julia-on-hpc-systems).

If you have [Homebrew](https://brew.sh) installed,
[open `Terminal.app`](https://support.apple.com/guide/terminal/open-or-quit-terminal-apd5265185d-f365-44cb-8b09-71a064a42125/mac)
and type

```bash
brew install julia
```

to install it as a [formula](https://docs.brew.sh/Formula-Cookbook).

If you are also using macOS and want to install it as a prebuilt binary app, type

```bash
brew install --cask julia
```

instead.

If you want to install multiple Julia versions in the same operating system,
a recommended way is to use a version manager such as
[`juliaup`](https://github.com/JuliaLang/juliaup).
First, [install `juliaup`](https://github.com/JuliaLang/juliaup#installation).
Then, run

```bash
juliaup add release
juliaup default release
```

to configure the `julia` command to start the latest stable version of
Julia (this is also the default value).

There is a [short video introduction to `juliaup`](https://youtu.be/14zfdbzq5BM)
made by its authors.

### Which version should I pick?

You can install the "Current stable release" or the "Long-term support (LTS)
release".

- The "Current stable release" is the latest release of Julia. It has access to
  newer features, and is likely faster.
- The "Long-term support release" is an older version of Julia that has
  continued to receive bug and security fixes. However, it may not have the
  latest features or performance improvements.

For most users, you should install the "Current stable release", and whenever
Julia releases a new version of the current stable release, you should update
your version of Julia. Note that any code you write on one version of the
current stable release will continue to work on all subsequent releases.

For users in restricted software environments (e.g., your enterprise IT controls
what software you can install), you may be better off installing the long-term
support release because you will not have to update Julia as frequently.

Versions higher than `v1.3`,
especially `v1.6`, are strongly recommended. This package may not work on `v1.0` and below.
Since the Julia team has set `v1.6` as the LTS release,
we will gradually drop support for versions below `v1.6`.

Julia and Julia packages support multiple operating systems and CPU architectures; check
[this table](https://julialang.org/downloads/#supported_platforms) to see if it can be
installed on your machine. For Mac computers with M-series processors, this package and its
dependencies may not work. Please install the Intel-compatible version of Julia (for macOS
x86-64) if any platform-related error occurs.

## Install ExpressBase

Now I am using [macOS](https://en.wikipedia.org/wiki/MacOS) as a standard
platform to explain the following steps:

1. Open `Terminal.app`, and type `julia` to start an interactive session (known as the
   [REPL](https://docs.julialang.org/en/v1/stdlib/REPL/)).

2. Run the following commands and wait for them to finish:

   ```julia-repl
   julia> using Pkg

   julia> Pkg.update()

   julia> Pkg.add("ExpressBase")
   ```

3. Run

   ```julia-repl
   julia> using ExpressBase
   ```

   and have fun!

4. While using, please keep this Julia session alive. Restarting might cost some time.

If you want to install the latest in-development (probably buggy)
version of ExpressBase, type

```@repl
using Pkg
Pkg.update()
pkg"add https://github.com/MineralsCloud/ExpressBase.jl"
```

in the second step above.

## Update ExpressBase

Please [watch](https://docs.github.com/en/account-and-profile/managing-subscriptions-and-notifications-on-github/setting-up-notifications/configuring-notifications#configuring-your-watch-settings-for-an-individual-repository)
our [GitHub repository](https://github.com/MineralsCloud/ExpressBase.jl)
for new releases.
Once we release a new version, you can update ExpressBase by typing

```@repl
using Pkg
Pkg.update("ExpressBase")
Pkg.gc()
```

in the Julia REPL.

## Uninstall and reinstall ExpressBase

Sometimes errors may occur if the package is not properly installed.
In this case, you may want to uninstall and reinstall the package. Here is how to do that:

1. To uninstall, in a Julia session, run

   ```julia-repl
   julia> using Pkg

   julia> Pkg.rm("ExpressBase")

   julia> Pkg.gc()
   ```

2. Press `ctrl+d` to quit the current session. Start a new Julia session and
   reinstall ExpressBase.
