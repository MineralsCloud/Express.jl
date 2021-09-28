<div align="center">
  <img src="docs/src/assets/logo.png" height="400"><br>
</div>

# Express

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://MineralsCloud.github.io/Express.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://MineralsCloud.github.io/Express.jl/dev)
[![Build Status](https://github.com/MineralsCloud/Express.jl/workflows/CI/badge.svg)](https://github.com/MineralsCloud/Express.jl/actions)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/MineralsCloud/Express.jl?svg=true)](https://ci.appveyor.com/project/singularitti/Express-jl)
[![Build Status](https://cloud.drone.io/api/badges/MineralsCloud/Express.jl/status.svg)](https://cloud.drone.io/MineralsCloud/Express.jl)
[![Build Status](https://api.cirrus-ci.com/github/MineralsCloud/Express.jl.svg)](https://cirrus-ci.com/github/MineralsCloud/Express.jl)
[![pipeline status](https://gitlab.com/singularitti/Express.jl/badges/master/pipeline.svg)](https://gitlab.com/singularitti/Express.jl/-/pipelines)
[![coverage report](https://gitlab.com/singularitti/Express.jl/badges/master/coverage.svg)](https://gitlab.com/singularitti/Express.jl/-/jobs)
[![Coverage](https://codecov.io/gh/MineralsCloud/Express.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/MineralsCloud/Express.jl)
[![Open in Visual Studio Code](https://open.vscode.dev/badges/open-in-vscode.svg)](https://open.vscode.dev/organization/repository)

Please cite this article as:
```bibtex
@misc{zhang2021textttexpress,
      title={$\texttt{express}$: extensible, high-level workflows for swifter $\textit{ab initio}$ materials modeling}, 
      author={Qi Zhang and Chaoxuan Gu and Jingyi Zhuang and Renata M. Wentzcovitch},
      year={2021},
      eprint={2109.11724},
      archivePrefix={arXiv},
      primaryClass={physics.comp-ph}
}
```

`Express` is an extensible, high-throughput, high-level workflow framework that aims to
automate *ab initio* calculations for the materials science community. `Express` is shipped
with well-tested workflow templates, including structure optimization, equation of state
(EOS) fitting, phonon spectrum (lattice dynamics) calculation, and thermodynamic property
calculation in the framework of the quasi-harmonic approximation (QHA). It is designed to be
highly modularized so that its components can be reused across various occasions, and
customized workflows can be built on top of that.

Please go to our website for the most updated version of
[code](https://github.com/MineralsCloud/Express.jl) and
[documentation](https://mineralscloud.github.io/Express.jl/dev/).

The components of `express`:

![components](docs/src/assets/components.png)

- [`Express.jl`](https://github.com/MineralsCloud/Express.jl) provides a high-level
  interface to all the workflows, including file reading and writing, job creation,
  submission, monitoring, result retrieving, and data analysis. To work with specific
  software, install the corresponding plugin, e.g., `QuantumESPRESSOExpress.jl` for Quantum
  ESPRESSO.

- [`ExpressCommands.jl`](https://github.com/MineralsCloud/ExpressCommands.jl) is a
  user-friendly command-line interface of `Express.jl` for non-developers. It installs an
  executable '`xps`' that can execute code from configuration files provided by users.

- [`EquationsOfStateOfSolids.jl`](https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl)
  fits energy (or pressure) vs. volume results to equations of state, etc. These features
  are repetitively used in the equation of state workflow.

- [`Crystallography.jl`](https://github.com/MineralsCloud/Crystallography.jl) calculates a
  crystal's primitive cell (or supercell) volume from lattice parameters, finds symmetry
  operations and generates high symmetry points in the Brillouin zone, etc.

- [`PyQHA.jl`](https://github.com/MineralsCloud/PyQHA.jl) is a `Julia` wrapper of the
  `Python` [`qha` package](https://github.com/MineralsCloud/qha), which can calculate
  several thermodynamic properties of both single- and multi-configuration crystalline
  materials in the framework of quasi-harmonic approximation (QHA). The `qha` code is the
  foundation of the QHA workflow.

- [`Geotherm.jl`](https://github.com/MineralsCloud/Geotherm.jl) is a `Julia` interpretation
  of the `Fortran` code we used in
  [this paper](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017GL073294), which
  calculates the isentropic temperature/pressure gradient (geotherm) using thermodynamic
  properties obtained with the QHA workflow.

- [`Pseudopotentials.jl`](https://github.com/MineralsCloud/Pseudopotentials.jl) presents a
  database for storing and querying pseudopotentials used in *ab initio* calculations.

- [`SimpleWorkflows.jl`](https://github.com/MineralsCloud/SimpleWorkflows.jl) is the
  skeleton of the workflow system, which defines building blocks, composition rules, and
  operation order of workflows.

The
[`QuantumESPRESSOExpress.jl`](https://github.com/MineralsCloud/QuantumESPRESSOExpress.jl) is
a special type of package called "plugin of `express`" for handling *ab initio* software
such as Quantum ESPRESSO. Other plugins for other software are possible. The dependencies of
`QuantumESPRESSOExpress.jl` are listed below.

- [`AbInitioSoftwareBase.jl`](https://github.com/MineralsCloud/AbInitioSoftwareBase.jl)
  provides a standard API for some popular *ab initio* software such as Quantum ESPRESSO.

- [`QuantumESPRESSOBase.jl`](https://github.com/MineralsCloud/QuantumESPRESSOBase.jl)
  declares basic data types and methods for manipulating crystal structures, generating
  input files for Quantum ESPRESSO, error checking before running, etc.

- [`QuantumESPRESSOParser.jl`](https://github.com/MineralsCloud/QuantumESPRESSOParser.jl)
  parses the input or output files of Quantum ESPRESSO to extract and analyze data.

- [`QuantumESPRESSOFormatter.jl`](https://github.com/MineralsCloud/QuantumESPRESSOFormatter.jl)
  formats the input files of Quantum ESPRESSO.

- [`QuantumESPRESSOCommands.jl`](https://github.com/MineralsCloud/QuantumESPRESSOCommands.jl)
  is a command-line interface that exports the commands Quantum ESPRESSO uses in a
  configurable way.

- [`QuantumESPRESSO.jl`](https://github.com/MineralsCloud/QuantumESPRESSO.jl) is simply a
  wrapper of the types, methods, and commands defined in `QuantumESPRESSOBase.jl`,
  `QuantumESPRESSOParser.jl`, `QuantumESPRESSOFormatter.jl`, and
  `QuantumESPRESSOCommands.jl` under a common namespace.
