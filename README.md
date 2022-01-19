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
