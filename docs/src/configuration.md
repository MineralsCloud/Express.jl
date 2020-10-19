# Configuration files

`Express` can be run from a configuration file, with some preset rules.
The following sections introduce how to write such configuration files.
By now, only
[`YAML`](https://docs.ansible.com/ansible/latest/reference_appendices/YAMLSyntax.html),
[`JSON`](https://restfulapi.net/json-syntax/), and
[`TOML`](https://toml.io/en/) formats are supported. Please refer their official
documentation for their syntax.
In the examples below, I will use `YAML` syntax for configuration files.

## Fitting equations of state

- The first key `workflow` means which workflow you want to apply. Currently, only `eos` and
  `eos+phonon` are available.
- The next key is the software being used. If using Quantum ESPRESSO, write `qe`.
  - `bin`: means the path of the binary executable to run actually calculation. It can be
    a full path if the location of the binary is not in users `$PATH` environment variable.
  - `manager` means the process manager to be used. By now, only set `local`.
  - `n` is the total number of processes available. If running on a node with 24
    processors, set it to `24`. If multiple nodes are being used, just add them up.
- `dir` is the working directory. When running, the directory will be changed to there. That
  is, if you are using relative paths for `template`, `bin`, etc., the root directory will
  be `dir`.
- `template` is the template file for scf and vc-relax calculations. It is basically an
  input file for Quantum ESPRESSO, VASP, etc.
- `pressures` are the pressures that on the compression curve, they are usually the desired
  pressures for further calculations. The unit of them is gigapascal. It is usually standard
  to have at least 6 pressures and at least 1 negative pressure.
- `trial_eos` is the starting equation of state for setting volumes for corresponding
  pressures.
  - `type` is the type of that equation of state. Available options are `m` (Murnaghan EOS),
    `bm2` - `bm4` (Birch--Murnaghan second to fourth order EOSs) and `v` (Vinet EOS).
  - `parameters` are the parameters of that equation of state. With the first parameter be
    zero-pressure volume (``V_0``), the second be zero-pressure bulk modulus (``B_0``),
    the third be zero-pressure bulk modulus derivative (``B_0'``).
  - `units` are the units of the corresponding parameters. Allowed values for ``V_0`` are
    `angstrom^3`, `bohr^3`, `nm^3`, `pm^3`, etc. Allowed values for ``B_0`` are
    `GPa`, `Pa`, `Mbar`, `kbar`, `eV/angstrom^3`, `eV/bohr^3`, `eV/nm^3`,
    `Ry/angstrom^3`, `Ry/bohr^3`, `hartree/angstrom^3`, etc. ``B_0'`` is a dimensionless
    number so its unit must be `"1"` (Note the quotation marks.).

```yaml
workflow: eos
qe:
  bin: pw.x
  manager: local
  n: 16
dir: examples/GaN
template: examples/GaN/template.in
pressures:
  - -5
  - 0
  - 5
  - 10
  - 15
  - 20
  - 25
  - 30
trial_eos:
  type: bm3
  parameters:
    - 317
    - 210
    - 4
  units:
    - bohr^3
    - GPa
    - "1"
```