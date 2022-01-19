# Configuration files

`Express` can be run from a configuration file, with some preset rules.
The following sections introduce how to write such configuration files.
By now, only
[YAML](https://docs.ansible.com/ansible/latest/reference_appendices/YAMLSyntax.html),
[JSON](https://restfulapi.net/json-syntax/), and
[TOML](https://toml.io/en/) formats are supported. Please refer their official
documentation for their syntax.
In the examples below, we will use YAML syntax for configuration files.
But for readability purposes, we suggest users use the TOML syntax.

## Fitting equations of state

The configuration file for the EOS workflow has the following syntax:

- `recipe`: A string that represents the type of the workflow. Allowed value is `eos`.
- `template`: The path to a template input file for a specific software. It should be on the
  same file system where `express` is deployed.
- `trial_eos`: The trial EOS contains initial values for input files generation and EOS fitting.
  - `type`: A string that represents the type of the EOS. Allowed values are `murnaghan`
    (Murnaghan), `bm2` (Birch--Murnaghan second order), `bm3`, `bm4`, `vinet` (Vinet), `pt2`
    (Poirier--Tarantola second order), `pt3`, and `pt4`.
  - `values`: A vector of strings that specifies each value of the EOS.
    The default order is ``V_0``, ``B_0``, ``B'_0``(, ``B''_0``, etc.). Units must be provided.
- `fixed`:
  - `pressures` or `volumes`: Whether to fix pressures of volumes.
    - `values`: Specify the pressures or volumes. It can be a vector of numbers, or a string
      with the syntax `start:step:stop` to form an arithmetic sequence where `start`,
      `stop`, and `step` are numbers indicating the start, the end, and the common
      difference of that sequence. See
      ["Creating arrays using range objects"](https://en.wikibooks.org/wiki/Introducing_Julia/Arrays_and_tuples#Creating_arrays_using_range_objects)
      for more information.
    - `unit`: The units of pressure or volume. The pressure and volume default units are
      `GPa` and `angstrom^3`. Allowed values for volumes are `nm^3`, `angstrom^3`, `bohr^3`,
      etc. Allowed values for pressures are `Pa`, `GPa`, `TPa`, ..., `bar`, `kbar`, ...,
      `atm`, and the combinations of `eV`, `Ry`, `hartree`, `J`, with any unit of volume.
- `files`:
  - `dirs`: It specifies the paths of output directories.
    - `root`: The path of the root directory of output files.
    - `pattern`: A string specifying the naming convention of the output directories. Its
      default value is `p=`. For example, if `fixed.pressures.values` is a vector of
      pressures `[10, 20, 30]` which represents the relaxations are done from ``10-30``GPa,
      then the generated inputs and outputs will be stored in directories `p=10`, `p=20` and
      `p=30`.
- `save`:
  - `status`: The path to a binary file that stores the status of the workflow.
  - `eos`: The path to a binary file that stores the fitted equations of state.
- `cli`:
  - `mpi`: The configurations of the MPI software.
  - `np`: An integer indicating the number of processors/cores/CPUs used.

The code block below shows a typical configuration file for an EOS workflow in the YAML syntax:

```yaml
recipe: eos
cli:
  mpi:
    np: 16
template: template.in
save:
  status: status.jls
fixed:
  pressures:
    unit: GPa
    values:
      - -5
      - -2
      - 0
      - 5
      - 10
      - 15
      - 17
      - 20
trial_eos:
  type: bm3
  values:
    - 300.44 bohr^3
    - 74.88 GPa
    - 4.82
```

The JSON and TOML equivalents of the above file are:

```json
{
  "recipe": "eos",
  "cli": {
    "mpi": {
      "np": 16
    }
  },
  "template": "template.in",
  "save": {
    "status": "status.jls"
  },
  "fixed": {
    "pressures": {
      "unit": "GPa",
      "values": [-5, -2, 0, 5, 10, 15, 17, 20]
    }
  },
  "trial_eos": {
    "type": "bm3",
    "values": ["300.44 bohr^3", "74.88 GPa", 4.82]
  }
}
```

```toml
recipe = 'eos'
template = 'template.in'
[fixed.pressures]
unit = 'GPa'
values = [-5, -2, 0, 5, 10, 15, 17, 20]

[trial_eos]
type = 'bm3'
values = ['300.44 bohr^3', '74.88 GPa', 4.82]
[cli.mpi]
np = 16

[save]
status = 'status.jls'
```

## Phonon density of states or phonon dispersion relation

The configuration file for the phonon workflow has the following syntax:

- `recipe`: A string that represents the type of the workflow. Allowed values are
  `phonon dispersion` (phonon dispersion along a q-path) and `vdos` (phonon density of states).
- `template`:
  - `scf`: The path to a template input file for an SCF calculation.
  - `dfpt`: The path to a template input file for a DFPT calculation.
  - `q2r`: The path to a template input file for a Fourier transform.
  - `disp`: The path to a template input file for a phonon dispersion/phonon density of states calculation.
- `fixed`:
  - `pressures` or `volumes`: Whether to fix pressures of volumes.
    - `values`: Specify the pressures or volumes. It can be a vector of numbers, or a string
      with the syntax `start:step:stop` to form an arithmetic sequence where `start`,
      `stop`, and `step` are numbers indicating the start, the end, and the common
      difference of that sequence. See
      ["Creating arrays using range objects"](https://en.wikibooks.org/wiki/Introducing_Julia/Arrays_and_tuples#Creating_arrays_using_range_objects)
      for more information.
    - `unit`: The units of pressure or volume. The pressure and volume default units are
      `GPa` and `angstrom^3`. Allowed values for volumes are `nm^3`, `angstrom^3`, `bohr^3`,
      etc. Allowed values for pressures are `Pa`, `GPa`, `TPa`, ..., `bar`, `kbar`, ...,
      `atm`, and the combinations of `eV`, `Ry`, `hartree`, `J`, with any unit of volume.
- `files`:
  - `dirs`: It specifies the paths of output directories.
    - `root`: The path of the root directory of output files.
    - `pattern`: A string specifying the naming convention of the output directories. Its
      default value is `p=`. For example, if `fixed.pressures.values` is a vector of
      pressures `[10, 20, 30]` which represents the relaxations are done from ``10-30``GPa,
      then the generated inputs and outputs will be stored in directories `p=10`, `p=20` and
      `p=30`.
- `save`:
  - `status`: The path to a binary file that stores the status of the workflow.
- `cli`:
  - `mpi`: The configurations of the MPI software.
  - `np`: An integer indicating the number of processors/cores/CPUs used.

The code block below shows a typical configuration file for a phonon workflow in the YAML syntax:

```yaml
recipe: vdos
cli:
  mpi:
    np: 16
template:
  scf: template.in
  dfpt: ph.in
  q2r: q2r.in
  disp: disp.in
save:
  status: status.jls
fixed:
  pressures:
    unit: GPa
    values:
      - -5
      - -2
      - 0
      - 5
      - 10
      - 15
      - 17
      - 20
```

The JSON and TOML equivalents of the above file are:

```json
{
  "recipe": "vdos",
  "cli": {
    "mpi": {
      "np": 16
    }
  },
  "template": {
    "scf": "template.in",
    "dfpt": "ph.in",
    "q2r": "q2r.in",
    "disp": "disp.in"
  },
  "save": {
    "status": "status.jls"
  },
  "fixed": {
    "pressures": {
      "unit": "GPa",
      "values": [-5, -2, 0, 5, 10, 15, 17, 20]
    }
  }
}
```

```toml
recipe = 'vdos'
[cli.mpi]
np = 16

[template]
dfpt = 'ph.in'
disp = 'disp.in'
q2r = 'q2r.in'
scf = 'template.in'
[fixed.pressures]
unit = 'GPa'
values = [-5, -2, 0, 5, 10, 15, 17, 20]

[save]
status = 'status.jls'
```
