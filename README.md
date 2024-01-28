<div align="center">
  <img src="https://raw.githubusercontent.com/MineralsCloud/Express.jl/master/docs/src/assets/logo.png" height="400"><br>
</div>

# Express: a high-level, extensible workflow framework for accelerating _ab initio_ calculations

|                                 **Documentation**                                  |                                                                                                 **Build Status**                                                                                                 |                                        **Others**                                         |
| :--------------------------------------------------------------------------------: | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: | :---------------------------------------------------------------------------------------: |
| [![Stable][docs-stable-img]][docs-stable-url] [![Dev][docs-dev-img]][docs-dev-url] | [![Build Status][gha-img]][gha-url] [![Build Status][appveyor-img]][appveyor-url] [![Build Status][cirrus-img]][cirrus-url] [![pipeline status][gitlab-img]][gitlab-url] [![Coverage][codecov-img]][codecov-url] | [![GitHub license][license-img]][license-url] [![Code Style: Blue][style-img]][style-url] |

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://MineralsCloud.github.io/Express.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://MineralsCloud.github.io/Express.jl/dev
[gha-img]: https://github.com/MineralsCloud/Express.jl/workflows/CI/badge.svg
[gha-url]: https://github.com/MineralsCloud/Express.jl/actions
[appveyor-img]: https://ci.appveyor.com/api/projects/status/github/MineralsCloud/Express.jl?svg=true
[appveyor-url]: https://ci.appveyor.com/project/singularitti/Express-jl
[cirrus-img]: https://api.cirrus-ci.com/github/MineralsCloud/Express.jl.svg
[cirrus-url]: https://cirrus-ci.com/github/MineralsCloud/Express.jl
[gitlab-img]: https://gitlab.com/singularitti/Express.jl/badges/main/pipeline.svg
[gitlab-url]: https://gitlab.com/singularitti/Express.jl/-/pipelines
[codecov-img]: https://codecov.io/gh/MineralsCloud/Express.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/MineralsCloud/Express.jl
[license-img]: https://img.shields.io/github/license/MineralsCloud/Express.jl
[license-url]: https://github.com/MineralsCloud/Express.jl/blob/main/LICENSE
[style-img]: https://img.shields.io/badge/code%20style-blue-4495d1.svg
[style-url]: https://github.com/invenia/BlueStyle

"Express" is a extensible, high-performance workflow framework designed to automate _ab initio_
calculations in the field of materials science. This framework includes rigorously tested
workflow templates for tasks such as structure optimization, equation of state (EOS)
fitting, phonon spectrum (lattice dynamics) computation, and thermodynamic property
calculation within the quasi-harmonic approximation (QHA) framework. Express is structured
with modularity in mind, enabling the reuse of its components across diverse applications
and facilitating the development of customized workflows.

The code is [hosted on GitHub](https://github.com/MineralsCloud/Express.jl),
with some continuous integration services to test its validity.

This repository is created and maintained by [@singularitti](https://github.com/singularitti).
You are very welcome to contribute.

Please [cite this package](https://doi.org/10.1016/j.cpc.2022.108515) as:

Q. Zhang, C. Gu, J. Zhuang et al., `express`: extensible, high-level workflows for swifter *ab initio* materials modeling, *Computer Physics Communications*, 108515, doi: https://doi.org/10.1016/j.cpc.2022.108515.

The BibTeX format is:

```bibtex
@article{ZHANG2022108515,
  title    = {express: extensible, high-level workflows for swifter ab initio materials modeling},
  journal  = {Computer Physics Communications},
  pages    = {108515},
  year     = {2022},
  issn     = {0010-4655},
  doi      = {https://doi.org/10.1016/j.cpc.2022.108515},
  url      = {https://www.sciencedirect.com/science/article/pii/S001046552200234X},
  author   = {Qi Zhang and Chaoxuan Gu and Jingyi Zhuang and Renata M. Wentzcovitch},
  keywords = {automation, workflow, high-level, high-throughput, data lineage}
}
```

We also have an [arXiv prepint](https://arxiv.org/abs/2109.11724).

## Installation

The package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
pkg> add Express
```

Or, equivalently, via the [`Pkg` API](https://pkgdocs.julialang.org/v1/getting-started/):

```julia
julia> import Pkg; Pkg.add("Express")
```

## Documentation

- [**STABLE**][docs-stable-url] — **documentation of the most recently tagged version.**
- [**DEV**][docs-dev-url] — _documentation of the in-development version._

Here is a video introduction to this code:

[![Watch the video](https://img.youtube.com/vi/N5_NUIaXnng/maxresdefault.jpg)](https://youtu.be/N5_NUIaXnng)

## Project status

The package is tested against, and being developed for, Julia `1.6` and above on Linux,
macOS, and Windows.

## Questions and contributions

You are welcome to post usage questions on [our discussion page][discussions-url].

Contributions are very welcome, as are feature requests and suggestions. Please open an
[issue][issues-url] if you encounter any problems. The [Contributing](@ref) page has
guidelines that should be followed when opening pull requests and contributing code.

## Stargazers over time

[![Stargazers over time](https://starchart.cc/MineralsCloud/Express.jl.svg?variant=adaptive)](https://starchart.cc/MineralsCloud/Express.jl)

[discussions-url]: https://github.com/MineralsCloud/Express.jl/discussions
[issues-url]: https://github.com/MineralsCloud/Express.jl/issues
