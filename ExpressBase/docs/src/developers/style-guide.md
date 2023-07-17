# [Style Guide](@id style)

This section describes the coding style rules that apply to our code and that
we recommend you to use it also.

In some cases, our style guide diverges from Julia's official
[Style Guide](https://docs.julialang.org/en/v1/manual/style-guide/) (Please read it!).
All such cases will be explicitly noted and justified.

Our style guide adopts many recommendations from the
[BlueStyle](https://github.com/invenia/BlueStyle).
Please read the [BlueStyle](https://github.com/invenia/BlueStyle)
before contributing to this package.
If not following, your pull requests may not be accepted.

!!! info
    The style guide is always a work in progress, and not all ExpressBase code
    follows the rules. When modifying ExpressBase, please fix the style violations
    of the surrounding code (i.e., leave the code tidier than when you
    started). If large changes are needed, consider separating them into
    another pull request.

## Formatting

### Run JuliaFormatter

ExpressBase uses [JuliaFormatter](https://github.com/domluna/JuliaFormatter.jl) as
an auto-formatting tool.

We use the options contained in [`.JuliaFormatter.toml`](https://github.com/MineralsCloud/ExpressBase.jl/blob/main/.JuliaFormatter.toml).

To format your code, `cd` to the ExpressBase directory, then run:

```@repl
using Pkg
Pkg.add("JuliaFormatter")
using JuliaFormatter: format
format("docs")
format("src")
format("test")
```

!!! info
    A continuous integration check verifies that all PRs made to ExpressBase have
    passed the formatter.

The following sections outline extra style guide points that are not fixed
automatically by JuliaFormatter.

### Use the Julia extension for Visual Studio Code

Please use [Visual Studio Code](https://code.visualstudio.com/) with the
[Julia extension](https://marketplace.visualstudio.com/items?itemName=julialang.language-julia)
to edit, format, and test your code.
We do not recommend using other editors to edit your code for the time being.

This extension already has [JuliaFormatter](https://github.com/domluna/JuliaFormatter.jl)
integrated. So to format your code, follow the steps listed
[here](https://www.julia-vscode.org/docs/stable/userguide/formatter/).
