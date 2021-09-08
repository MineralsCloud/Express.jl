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

- Next, specify how many computing cores you are using in total. If running on multiple nodes,
  write the summation of cores on these nodes. It will be better if `np` is an integer
  multiple of the number of pressures.
- The next key `bin` is the software being used. If using Quantum ESPRESSO, write `qe`. The
  values of `qe` is the path to the binaries that actually run the calculations. If they are
  already in the `PATH` environment variable, write the names of the executables.
- `templates` is a vector of the template files for scf and vc-relax calculations. If it has
  only one value, all the pressures share the same template. If it has more than one value,
  the number of files must be equal to the number of pressures. That is, each pressure has
  its own template.
- `pressures` are the pressures that on the compression curve, they are usually the desired
  pressures for further calculations. The unit of them is by default `GPa`. It is usually
  standard to have at least 6 pressures and at least 1 negative pressure.
- `trial_eos` is the starting equation of state for setting volumes for corresponding
  pressures.
  - `name` is the name of that equation of state. Available options are `m` (Murnaghan EOS),
    `bm2` - `bm4` (Birch--Murnaghan second to fourth order EOSs) and `v` (Vinet EOS).
  - `parameters` are the parameters of that equation of state. With the first parameter be
    zero-pressure volume (``V_0``), the second be zero-pressure bulk modulus (``B_0``), the
    third be zero-pressure bulk modulus derivative (``B_0'``). Each value has an associate
    unit. Allowed units for ``V_0`` are `angstrom^3`, `bohr^3`, `nm^3`, `pm^3`, etc. Allowed
    values for ``B_0`` are `GPa`, `Pa`, `Mbar`, `kbar`, `eV/angstrom^3`, `eV/bohr^3`,
    `eV/nm^3`, `Ry/angstrom^3`, `Ry/bohr^3`, `hartree/angstrom^3`, etc. ``B_0'`` is a
    dimensionless number so its unit **must be** `1`.
- `use_shell`: Whether create shell files to run external software to do computations.
  Usually this is preferred when `express` is run non-interactively. If run in interactive
  mode, you may want to set it to `false`.

```yaml
workflow: eos
np: 24
bin:
  qe:
    - pw.x
templates:
  - examples/Ge/template.in
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
  name: bm3
  parameters:
    - - 300.44
      - bohr^3
    - - 74.88
      - GPa
    - - 4.82
      - 1
    - - -612.43149513
      - Ry
use_shell: true
```

## Phonon density of states or phonon dispersion relation

- The first key `workflow` means which workflow you want to apply. For Phonon density of states,
  write `vdos`. For phonon dispersion relation, write `phonon dispersion`.
- Next, specify how many computing cores you are using in total. If running on multiple nodes,
  write the summation of cores on these nodes. It will be better if `np` is an integer
  multiple of the number of pressures.
- The next key `bin` is the software being used. If using Quantum ESPRESSO, write `qe`. The
  values of `qe` is the path to the binaries that actually run the calculations. If they are
  already in the `PATH` environment variable, write the names of the executables. They
  should be written sequentially, with each one corresponds to one calculation.
- `templates` is a vector of the template files for the calculations. If it has
  only one value, all the pressures share the same template. If it has more than one value,
  the number of files must be equal to the number of pressures. That is, each pressure has
  its own template. Note each calculation has a vector of template files.
- `pressures` are the pressures that on the compression curve, they are usually the desired
  pressures for further calculations. The unit of them is by default `GPa`. It is usually
  standard to have at least 6 pressures and at least 1 negative pressure.
- `use_shell`: Whether create shell files to run external software to do computations.
  Usually this is preferred when `express` is run non-interactively. If run in interactive
  mode, you may want to set it to `false`.

```yaml
workflow: vdos
np: 24
bin:
  qe:
    - pw.x
    - ph.x
    - q2r.x
    - matdyn.x
templates:
  - - examples/Ge/template.in
  - - examples/Ge/ph.in
  - - examples/Ge/q2r.in
  - - examples/Ge/disp.in
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
use_shell: true
```
