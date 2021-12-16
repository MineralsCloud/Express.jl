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

Please view the [arXiv prepint here](https://arxiv.org/abs/2109.11724) and cite this package as:
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

- [`PyQHA.jl`](https://github.com/MineralsCloud/PyQHA.jl) is a Julia wrapper of the
  Python [`qha` package](https://github.com/MineralsCloud/qha), which can calculate
  several thermodynamic properties of both single- and multi-configuration crystalline
  materials in the framework of quasi-harmonic approximation (QHA). The `qha` code is the
  foundation of the QHA workflow.

- [`Geotherm.jl`](https://github.com/MineralsCloud/Geotherm.jl) is a Julia interpretation
  of the Fortran code we used in
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

## Some questions about `Express`

### What is the difference between `Express` and `express`?

`express` is the workflow framework's name, while `Express` is the short form of
`Express.jl`, the Julia implementation of `express`. We do not want the project's name
linked to a specific programming language: who says we cannot have a Python version of
`express` in the future?

### Why do you use Julia, why not use...?

The same as the Julia's ["we are greedy" speech](https://julialang.org/blog/2012/02/why-we-created-julia/):

> We are greedy: we want more.
>
> We want a language that's open source, with a liberal license. We want the speed of C with
> the dynamism of Ruby. We want a language that's homoiconic, with true macros like Lisp,
> but with obvious, familiar mathematical notation like Matlab. We want something as usable
> for general programming as Python, as easy for statistics as R, as natural for string
> processing as Perl, as powerful for linear algebra as Matlab, as good at gluing programs
> together as the shell. Something that is dirt simple to learn, yet keeps the most serious
> hackers happy. We want it interactive and we want it compiled.

To summarize,
- we want the code as easy to read/write as Python, but as fast as C;
- we want the code as extensible as possible;
- we want the flourishing ecosystem of Julia, dedicated to and developed by the scientific community itself.

### Why do you create `express`, given that [`AiiDA`](https://www.aiida.net/), [`ASE`](https://gitlab.com/ase/ase), [`atomate`](https://atomate.org/), etc., are already there?

Good question. The short answer is: we are never satisfied.

Those projects are of very high quality and are de facto standards of the materials
simulation community. They have much larger teams and longer history than us, so
it is unrealistic to say we are better than them in every aspect in the first few releases
of `express`. However, it does not mean we cannot have our features or advantages.
Our code has higher modularity, extensibility, and readability.

As explained in [the paper](https://arxiv.org/abs/2109.11724), `express` is a highly
modularized collection of packages. Each of them provides a succinct but complete
set of functionalities that is repeatedly used in materials modeling. Even without
other packages, each of them can still solve its dedicated questions. We faithfully follow
the Unix philosophy

> Write programs that do one thing and do it well.
> Write programs to work together.

when developing `express`. For example,
[`Pseudopotentials.jl`](https://github.com/MineralsCloud/Pseudopotentials.jl),
[`Crystallography.jl`](https://github.com/MineralsCloud/Crystallography.jl),
[`EquationsOfStateOfSolids.jl`](https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl),
[`QuantumESPRESSO.jl`](https://github.com/MineralsCloud/QuantumESPRESSO.jl), etc.,
are distinct packages. They can work together to form the equation of state workflow, and
they can also be used separately if users want to use specific features, e.g., managing
pseudopotentials, calculating structural symmetry, fitting existing data.
In that case, only related code will be installed, which saves disk space and time. However,
most code in, for instance, [`aiida-core`](https://github.com/aiidateam/aiida-core), is
dealing with servers, network, databases, web interfaces, YAML, etc, which makes it a
gigantic project that is hard to read, trace, and debug. And everything in `AiiDA` is bond
to its workflow system. Its data structure may not be compatible with other people's code.
It is difficult for users who are unfamiliar with their code structure to pick a subset of
features they need if they want to customize. In `express`, everything is very loosely
coupled to each other, no workflow needs to be composed if you just want to do a simple
thing (In fact, each workflow in `express` is just a collection of predefined wrappers of
functions that its dependencies have already provided. The core code of `Express.jl` is very
small.). Users can cooperate with others' packages with little effort thanks to the
composability that Julia enables. For example, in `express`, our users never need to convert
units. With the help of Julia's [`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl)
package, users can always throw literal units they obtained from experiments or papers
directly to `express` and do not need to consider whether they are compatible or not. This
is helpful when fitting equations of state, modifying crystal structures, etc. From
developers' side, we seldom need to adapt for others' code. Theirs and ours just "magically"
work together. You usually do not have such flexibility in other languages. Besides, even
the core calculation part of `AiiDA` and `atomate`, they use
[`pymatgen`](https://pymatgen.org/), another phenomenal code in the materials science
community. However, it is still an enormous code that integrates too many things together,
which limits its extensibility.

Due to the modularity we mentioned above, the `express` project is separated into many
packages. Each package has its own documentation and tests, and is released individually to
avoid updating the whole codebase whenever there is a bug fix or a feature enhancement. They
also have separated pull request pages for developers and skilled users to discuss and
collaborate. Julia's semantic versioning system manages their compatibility, i.e.,
compatible packages are downloaded automatically, and no human intervention is needed.

In addition, some of the workflows we shipped in `express` is uniquely developed by us. See
the introduction of the [`qha` package](https://doi.org/10.1016/j.cpc.2018.11.003), which
can calculate quasiharmonic free energy for multi-configuration systems. We also have some
workflows that will be integrated into `express` in the near future, including but not
limited to, the `pgm` (phonon gas model) workflow that was used to calculate thermodynamic
properties of
[ε-Fe with thermal electronic excitation effects on vibrational spectra](https://doi.org/10.1103/PhysRevB.103.144102),
and the [thermoelasticity workflow](https://doi.org/10.1016/j.cpc.2021.108067) based on the
[Wu–Wentzcovitch semi-analytical method (SAM)](https://doi.org/10.1103/PhysRevB.83.184115).
