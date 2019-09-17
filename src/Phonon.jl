"""
# module Phonon



# Examples

```jldoctest
julia>
```
"""
module Phonon

using Kaleido: @batchlens
using QuantumESPRESSOBase: to_qe
using QuantumESPRESSOBase.Inputs.PWscf: PWscfInput
using QuantumESPRESSOParsers.OutputParsers.PWscf: read_cell_parameters, read_atomic_positions
using Setfield: get, set, @lens

import ..Step
using ..SelfConsistentField: write_metadata

export update_structure, relay, prepare

"""
    update_structure(output::AbstractString, template::PWscfInput)

Read structure information from `output`, and update related fields of `template`.
"""
function update_structure(output::AbstractString, template::PWscfInput)
    cell_parameters = read_cell_parameters(output)
    atomic_positions = read_atomic_positions(output)  # TODO: Implement `read_atomic_positions`
    lenses = @batchlens(begin
        _.system.celldm ∘ _[$1]
        _.atomic_positions
        _.cell_parameters
    end)
    return set(template, lenses, (1, atomic_positions, cell_parameters))
end # function update_structure

# This is a helper function and should not be exported.
function _preset(template::PWscfInput)
    lenses = @batchlens(begin
        _.control.calculation  # Get the `template`'s `control.calculation` value
        _.control.verbosity    # Get the `template`'s `control.verbosity` value
        _.control.tstress      # Get the `template`'s `control.tstress` value
        _.control.tprnfor      # Get the `template`'s `control.tprnfor` value
    end)
    # Set the `template`'s values with...
    template = set(template, lenses, ("scf", "high", true, true))
    return isnothing(template.cell_parameters) ? autofill_cell_parameters(template) : template
end # function _preset
function _preset(template::PhononInput)
    lenses = @batchlens(begin
        _.phonon.verbosity  # Get the `template`'s `phonon.calculation` value
    end)
    # Set the `template`'s values with...
    template = set(template, lenses, ("high",))
end # function _preset

"""
    relay(from::PWscfInput, to::PhononInput)

Relay shared information from a `PWscfInput` to a `PhononInput`.

A `PWscfInput` before a `PhononInput` has the information of `outdir` and `prefix`. They must keep the same in a
phonon calculation.
"""
function relay(from::PWscfInput, to::PhononInput)
    pwlenses = @batchlens(begin
        _.control.outdir
        _.control.prefix
    end)
    phlenses = @batchlens(begin
        _.phonon.outdir
        _.phonon.prefix
    end)
    return set(from, phlenses, get(to, pwlenses))
end # function relay
function relay(from::PhononInput, to::Q2RInput)
    fildyn = @lens _.fildyn
    return set(to, @lens _.q2r ∘ fildyn, get(@lens _.phonon ∘ fildyn))
end # function relay
function relay(from::Q2RInput, to::MatdynInput)
    q2r_lenses = @batchlens(begin
        _.q2r.fildyn
        _.q2r.flfrc
        _.q2r.loto_2d
        _.q2r.zasr
    end)
    matdyn_lenses = @batchlens(begin
        _.matdyn.fildyn
        _.matdyn.flfrc
        _.matdyn.loto_2d
        _.matdyn.asr
    end)
    return set(to, matdyn_lenses, get(from, q2r_lenses))
end # function relay
function relay(from::PhononInput, to::DynmatInput)
    ph_lenses = @batchlens(begin
        _.phonon.asr
        _.phonon.fildyn
        _.phonon.amass
    end)
    dynmat_lenses = @batchlens(begin
        _.phonon.asr
        _.phonon.fildyn
        _.phonon.amass
    end)
    return set(to, dynmat_lenses, get(from, ph_lenses))
end # function relay

"""
    prepare(step::Step{1}, inputs, outputs, template, metadatafiles[, verbose::Bool = false])

Prepare input files of the first step of a phonon calculation.

# Arguments
- `step::Step{1}`: Denotes the first step of a phonon calculation.
- `inputs::AbstractVector{<:AbstractString}`: The input files of Quantum ESPRESSO of a phonon calculation. If not exist,
   they will be automatically created.
- `outputs::AbstractVector{<:AbstractString}`: The output files of Quantum ESPRESSO of a vc-relax calculation.
- `template::PWscfInput`:
- `metadatafiles::AbstractVector{<:AbstractString}`:
- `verbose::Bool = false`: control the format of input files, verbose or not.
"""
function prepare(
    step::Step{1},
    inputs::AbstractVector{<:AbstractString},
    outputs::AbstractVector{<:AbstractString},
    template::PWscfInput,
    metadatafiles::AbstractVector{<:AbstractString} = map(x -> splitext(x)[1] * ".json", inputs),
    verbose::Bool = false
)
    # Check parameters
    @assert(
        length(inputs) == length(outputs) == length(metadatafiles),
        "The inputs, outputs and the metadata files must be the same length!"
    )
    template = _preset(template)
    for (input, output, structure, metadatafile) in zip(inputs, outputs, metadatafiles)
        # Get a new `object` from the `template`, with its `alat` and `pressure` changed
        object = update_structure(output, template)
        write(input, to_qe(object, verbose = verbose))  # Write the `object` to a Quantum ESPRESSO input file
        write_metadata(metadatafile, input, template)
    end
end # function prepare
function prepare(
    step::Step{2},
    phonon_inputs::AbstractVector{<:AbstractString},
    pwscf_inputs::AbstractVector{<:AbstractString},
    template::PhononInput,
    verbose::Bool = false
)
    # Check parameters
    @assert(
        length(phonon_inputs) == length(pwscf_inputs),
        "The PWscf and the PHonon inputs files must have the same length!"
    )
    template = _preset(template)
    for (phonon_input, pwscf_input) in zip(phonon_inputs, pwscf_inputs)
        open(pwscf_input, "r") do io
            object = parse(PWscfInput, read(io, String))
        end
        template = relay(template, object)
        write(phonon_input, template)
    end
end # function prepare

end
