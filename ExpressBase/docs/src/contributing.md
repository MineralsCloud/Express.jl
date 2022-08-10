# [Contributing](@id contributing)

```@contents
Pages = ["contributing.md"]
```

## Download the project

Similar to [installation](@ref), run

```@repl
using Pkg
Pkg.update()
pkg"dev ExpressBase"
```

in the REPL.

Then the package will be cloned to your local machine at a path. On macOS, by default is
located at `~/.julia/dev/ExpressBase` unless you modify the `JULIA_DEPOT_PATH`
environment variable. (See [Julia's official documentation](http://docs.julialang.org/en/v1/manual/environment-variables/#JULIA_DEPOT_PATH-1)
on how to do this.) In the following text, we will call it `PKGROOT`.

## [Instantiate the project](@id instantiating)

Go to `PKGROOT`, start a new Julia session and run

```@repl
using Pkg
Pkg.instantiate()
```

## How to build docs

Usually, the up-to-state doc is available in
[here](https://MineralsCloud.github.io/ExpressBase.jl/dev), but there are cases
where users need to build the doc themselves.

After [instantiating](@ref) the project, go to `PKGROOT`, run

```bash
julia --color=yes docs/make.jl
```

in your terminal. In a while a folder `PKGROOT/docs/build` will appear. Open
`PKGROOT/docs/build/index.html` with your favorite browser and have fun!

## How to run tests

After [instantiating](@ref) the project, go to `PKGROOT`, run

```bash
julia --color=yes test/runtests.jl
```

in your terminal.

## Style Guide

Follow the style of the surrounding text when making changes. When adding new features
please try to stick to the following points whenever applicable.

### Julia

- 4-space indentation;
- modules spanning entire files should not be indented, but modules that have surrounding code should;
- do not manually align syntax such as `=` or `::` over adjacent lines;
- use `function ... end` when a method definition contains more than one top-level expression;
- related short-form method definitions don't need a new line between them;
- unrelated or long-form method definitions must have a blank line separating each one;
- surround all binary operators with whitespace except for `::`, `^`, and `:`;
- files containing a single `module ... end` must be named after the module;
- method arguments should be ordered based on the amount of usage within the method body;
- methods extended from other modules must follow their inherited argument order, not the above rule;
- explicit `return` should be preferred except in short-form method definitions;
- avoid dense expressions where possible e.g. prefer nested `if`s over complex nested `?`s;
- include a trailing `,` in vectors, tuples, or method calls that span several lines;
- do not use multiline comments (`#=` and `=#`);
- wrap long lines as near to 92 characters as possible, this includes docstrings;
- follow the standard naming conventions used in `Base`.

### Markdown

- Use unbalanced `#` headers, i.e. no `#` on the right-hand side of the header text;
- include a single blank line between top-level blocks;
- do _not_ hard wrap lines;
- use emphasis (`*`) and bold (`**`) sparingly;
- always use fenced code blocks instead of indented blocks;
- follow the conventions outlined in the Julia documentation page on documentation.
