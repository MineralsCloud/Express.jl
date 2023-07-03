# How to run `Express` from command line

To carry out a computation, first, users need to prepare the input files introduced in
section "[Configuration files](@ref)". Then if [the `xps` command](@ref cli) is installed,
run in terminal:

```shell
xps run eos.toml  # equation of state workflow
xps run vdos.toml  # phonon workflow
xps run qha.toml  # QHA workflow
```

As explained in "[Configuration files](@ref)", the file name (`eos.toml`, etc.)
of the configuration file does not matter. It is the `recipe` field that
determines the type of workflow.

It is often suggested to put the above commands in a script file with a header
and submit it to a workload manager like [Slurm](https://www.schedmd.com/).
However, users can also request an interactive session and run code Julia REPL.

```julia-repl
julia> using Express.EquationOfStateWorkflow.Recipes

julia> using QuantumESPRESSOExpress

julia> workflow = buildworkflow("eos.toml");

julia> run!(workflow);
```

For phonon and the QHA workflow, run `using Express.PhononWorkflow.Recipes`
and `using Express.QuasiHarmonicApproxWorkflow.Recipes` in the first step.

In an EOS workflow, the final results include outputs returned by Quantum ESPRESSO, a list
of raw data (volume-energy pairs), and a fitted EOS. If something goes wrong, the workflow
might be terminated. The state of the workflow (i.e., the status of each job) and errors
will be saved in a file for debugging. Once the bug is fixed,
run `xps run <path-to-config-file>` again, and `express` will retry the failed jobs.
To print either input or output data in a formatted, readable form, run

```shell
xps print <file-name>
```

where the allowed extensions of `<file-name>` are `.jls`, `.json`, `.yaml` or `.yml`, and
`.toml`. The last four extensions correspond to three human-readable data-serialization file
formats, i.e., [JSON](https://www.json.org/json-en.html), [YAML](https://yaml.org/), and
[TOML](https://toml.io/en/), while `.jls` is a binary serialization format only recognizable
to Julia. For example,

```shell
xps print eos.toml
xps print raw.json
xps print eos.jls
```

`express` can also plot some data, such as the fitted EOS applied to a certain range of
volumes along with the raw data. The corresponding command is

```shell
xps plot <file-name>
```

where `<file-name>` refers to the EOS binary file with extension `.jls`.

These are the three
most important commands of `express`. These catchy commands cover all the functionalities we
have promised, including but not limited to unit conversion, pseudopotential downloading,
input validation and generation, calculation monitoring, task distribution, gathering and
analysis of results, error handling, logging, and visualization. We hope they can facilitate
tedious work as much as possible.
