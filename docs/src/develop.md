# How to develop this package

## Download the project

First,
[add an SSH key to your GitHub account](https://docs.github.com/en/github/authenticating-to-github/adding-a-new-ssh-key-to-your-github-account).
Request or wait for the administrator of this repo to give you an invitation for
triage.

Then, similar to section "[Installation](@ref)", run

```@repl
using Pkg; Pkg.update()
pkg"dev Express"
```

Then the package will be cloned to your local machine at a path. On Unix-like systems, by
default is located at `~/.julia/dev/Express` unless you modify the
`JULIA_DEPOT_PATH` environment variable (See
[Julia's official documentation](https://docs.julialang.org/en/v1/manual/environment-variables/#JULIA_DEPOT_PATH)).
In the following text, we will call it `PKGROOT`.

## [Instantiate the project](@id instantiating)

Go to `PKGROOT`, start a new Julia session and run

```@repl
using Pkg; Pkg.instantiate()
```

## How to build docs

Usually, the up-to-state doc is available in
[here](https://MineralsCloud.github.io/Express.jl/dev), but there
are cases where users need to build the doc themselves.

After [instantiating](@ref) the project, go to `PKGROOT`, run (without the `$`
prompt)

```bash
$ julia --color=yes docs/make.jl
```

in your terminal. In a while a folder `PKGROOT/docs/build` will appear. Open
`PKGROOT/docs/build/index.html` with your favorite browser and have fun!

## How to run tests

After [instantiating](@ref) the project, go to `PKGROOT`, run (without the `$`
prompt)

```bash
$ julia --color=yes test/runtests.jl
```

in your terminal.
