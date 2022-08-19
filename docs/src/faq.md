# FAQ

## What is the difference between `Express` and `express`?

`express` is the workflow framework's name, while `Express` is the short form of
`Express.jl`, the Julia implementation of `express`. We do not want the project's name
linked to a specific programming language: who says we cannot have a Python version of
`express` in the future?

## Why do you create `express`?

Those projects are of very high quality and are de facto standards of the materials
simulation community. They have much larger teams and longer history than us, so
it is unrealistic to say we are better than them in every aspect in the first few releases
of `express`. However, it does not mean we cannot have our features or advantages.
Our code is better at modularity, extensibility, and readability.

First, some packages only implement workflows for specific software. By far, Quantum
ESPRESSO is only supported by a few packages. Considering the number of users in the Quantum ESPRESSO
community (1000-2000 citations per year), there is a great need for an advanced and eclectic
workflow ecosystem.

Second, most of their effort is put into dealing with servers, networks, databases, web
interfaces, file formats, etc., while the code's core logic takes up only a small part,
leaving gigantic packages that are hard to understand and integrate into users' code. On the
contrary, `express` is a highly modularized collection of
packages, with each of them providing a succinct, almost independent, complete set of
functionalities frequently used in materials modeling. Inside each package, every type and
function are also loosely coupled. Users can pick up pieces of our code and incorporate them
into theirs effortlessly.
It is also possible for users to build customized workflows. In
fact, each workflow in `Express.jl` is just a collection of predefined wrappers of functions
provided by its dependencies. Part of the `express` project can be installed and used if
necessary.
For example,
[`Pseudopotentials.jl`](https://github.com/MineralsCloud/Pseudopotentials.jl),
[`Crystallography.jl`](https://github.com/MineralsCloud/Crystallography.jl),
[`EquationsOfStateOfSolids.jl`](https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl), and
[`QuantumESPRESSO.jl`](https://github.com/MineralsCloud/QuantumESPRESSO.jl),
are distinct packages. They can work together to form the equation of state workflow.
They can also work separately, such as managing
pseudopotentials, calculating structural symmetry, fitting existing data.

In addition, many packages mentioned above are implemented in Python, a convenient language
when building a prototype project but not the most convenient one when developing a large
project due to readability and performance issues (For this reason, Dropbox, Inc. is writing
a static type checker [`mypy`](http://mypy-lang.org/) to make large Python programs easier
to read, and people are working on a faster Python implementation called
[`PyPy`](https://www.pypy.org/).). However, Julia has built-in support for type annotations,
making it extremely readable and performant. Namely, it is easy to write generic Julia code
without losing performance. Another benefit of using Julia is good compatibility between
codebases. Because of its
[multimethods feature](https://docs.julialang.org/en/v1/manual/methods/),
different Julia packages usually just "magically" work together.
For example, our users never need to convert units.
With the help of [`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl),
we can always write literal units obtained from experiments
and do not need to consider whether they are the same or not. This
is helpful when fitting equations of state, modifying crystal structures, etc.
We usually do not have such flexibility in other languages,
including Python, where sometimes wrapper code is needed. From the developers' perspectives,
we seldom need to adapt to others' code. From the users' perspectives, they can use
customized data structures without worrying about code incompatibility.

At last, some of the workflows we shipped in `express` are uniquely developed by us. See
the introduction of the [`qha` package](https://doi.org/10.1016/j.cpc.2018.11.003), which
can calculate quasiharmonic free energy for multi-configuration systems. We also have some
workflows that will be integrated into `express` shortly, including but not
limited to, the phonon gas model workflow that was used to calculate thermodynamic
properties of
[ε-Fe with thermal electronic excitation effects on vibrational spectra](https://doi.org/10.1103/PhysRevB.103.144102),
and the [thermoelasticity workflow](https://doi.org/10.1016/j.cpc.2021.108067) based on the
[Wu–Wentzcovitch semi-analytical method (SAM)](https://doi.org/10.1103/PhysRevB.83.184115).
