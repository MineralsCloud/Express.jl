# Some other questions about `Express`

## What is the difference between `Express` and `express`?

`express` is the workflow framework's name, while `Express` is the short form of
`Express.jl`, the Julia implementation of `express`. We do not want the project's name
linked to a specific programming language: who says we cannot have a Python version of
`express` in the future?

## Why do you create `express`, given that [`AiiDA`](https://www.aiida.net/), [`ASE`](https://gitlab.com/ase/ase), [`atomate`](https://atomate.org/), etc., are already there?

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
are distinct packages. They can work together to form the equation of state workflow, or
they can be used separately if users want to use specific features, e.g., managing
pseudopotentials, calculating structural symmetry, fitting existing data.
In that case, only related code will be installed, which saves disk space and time. However,
most code in, for instance, [`aiida-core`](https://github.com/aiidateam/aiida-core),
deals with servers, networks, databases, web interfaces, YAML, etc., making it a
gigantic project that is hard to read, trace, and debug. And everything in `AiiDA` is bond
to its workflow system. Its data structure may not be compatible with other people's code.
It is difficult for users unfamiliar with their code structure to pick a subset of
features they need if they want to customize. In `express`, everything is very loosely
coupled to each other, and no workflow needs to be composed if you just want to do a simple
thing (In fact, each workflow in `express` is just a collection of predefined wrappers of
functions that its dependencies have already provided. The core code of `Express.jl` is very
small.). Users can cooperate with others' packages with little effort, thanks to the
composability that Julia enables. For example, in `express`, our users never need to convert
units. With the help of Julia's [`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl)
package, users can always throw literal units they obtained from experiments or papers
directly to `express` and do not need to consider whether they are compatible or not. This
is helpful when fitting equations of state, modifying crystal structures, etc. From the
developers' side, we seldom need to adapt to others' code. Theirs and ours just "magically"
work together. You usually do not have such flexibility in other languages. Besides, even
the core calculation part of `AiiDA` and `atomate`,
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
