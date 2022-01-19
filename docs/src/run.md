# How to run `Express` from command line

To carry out a computation, first, users need to prepare the input files introduced in
section "[Configuration files](@ref)". Then if [the `xps` command](@ref cli) is installed,
run in terminal:

```shell
$ xps run eos.toml  # equation of state workflow
$ xps run vdos.toml  # phonon workflow
$ xps run qha.toml  # QHA workflow
```

As explained in "[Configuration files](@ref)", the file name (`eos.toml`, etc.)
of the configuration file does not matter. It is the `recipe` field that
determines the type of the workflow.
