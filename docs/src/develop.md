# How to develop this package

## Download the project

First,
[add an SSH key to your GitHub account](https://docs.github.com/en/github/authenticating-to-github/adding-a-new-ssh-key-to-your-github-account).
Request or wait for the administrator of this repo to give you an invitation for triage.

Then, similar to section "[Installation](@ref)", run

```julia
julia> using Pkg; Pkg.update()

julia> pkg"add git@github.com:MineralsCloud/AbInitioSoftwareBase.jl.git"

julia> pkg"add git@github.com:MineralsCloud/Pseudopotentials.jl.git"

julia> pkg"add git@github.com:MineralsCloud/Crystallography.jl.git"

julia> pkg"add git@github.com:singularitti/PyFortran90Namelists.jl.git"

julia> pkg"add git@github.com:MineralsCloud/QuantumESPRESSOBase.jl.git"

julia> pkg"add git@github.com:MineralsCloud/QuantumESPRESSOParsers.jl.git"

julia> pkg"add git@github.com:MineralsCloud/QuantumESPRESSO.jl.git"

julia> pkg"dev git@github.com:MineralsCloud/Express.jl.git"
```

Then the package will be cloned to your local machine at a path. On macOS, by default is
located at `~/.julia/dev/Express` unless you modify the `JULIA_DEPOT_PATH`
environment variable. (See
[Julia's official documentation](http://docs.julialang.org/en/v1/manual/environment-variables/#JULIA_DEPOT_PATH-1)
on how to do this.) In the following text, we will call it `PKGROOT`.

## [Instantiate the project](@id instantiating)

Go to `PKGROOT`, start a new Julia session and run

```julia
julia> using Pkg; Pkg.instantiate()
```

## How to build docs

Usually, the up-to-state doc is available in
[here](https://MineralsCloud.github.io/AbInitioSoftwareBase.jl/dev), but there are cases
where users need to build the doc themselves.

After [instantiating](@ref) the project, go to `PKGROOT`, run (without the `$` prompt)

```bash
$ julia --color=yes docs/make.jl
```

in your terminal. In a while a folder `PKGROOT/docs/build` will appear. Open
`PKGROOT/docs/build/index.html` with your favorite browser and have fun!

## How to run tests

After [instantiating](@ref) the project, go to `PKGROOT`, run (without the `$` prompt)

```bash
$ julia --color=yes test/runtests.jl
```

in your terminal.
