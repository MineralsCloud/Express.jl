"""
# module Phonon



# Examples

```jldoctest
julia>
```
"""
module Phonon

using Kaleido: @batchlens
using QuantumESPRESSO: to_qe
using QuantumESPRESSO.Cards: option_convert
using QuantumESPRESSO.Cards.PWscf: AtomicPositionsCard, CellParametersCard
using QuantumESPRESSO.Cards.PHonon: SpecialQPoint, QPointsSpecsCard
using QuantumESPRESSO.Inputs.PWscf: PWInput
using QuantumESPRESSO.Inputs.PHonon: PhInput, Q2rInput, MatdynInput, DynmatInput
using QuantumESPRESSO.Outputs.PWscf: parsefinal
using QuantumESPRESSOBase.Setters: CellParametersSetter, batchset
using Setfield: get, set, @lens, @set!

import ..Step
using Express.BandStructure: generate_path

export update_structure, relay, preprocess

"""
    update_structure(output::AbstractString, template::PWInput)

Read structure information from `output`, and update related fields of `template`.
"""
function update_structure(output::AbstractString, template::PWInput)
    open(output, "r") do io
        str = read(io, String)
        # The result of `parsefinal` must be a `CellParametersCard` with `"bohr"` or `"angstrom"` option.
        # Convert it to "bohr" by default
        cell_parameters = option_convert("bohr", parsefinal(CellParametersCard, str))
        atomic_positions = parsefinal(AtomicPositionsCard, str)
        lenses = @batchlens(begin
            _.system.ibrav
            _.system.celldm ∘ _[$1]
            _.atomic_positions
            _.cell_parameters
        end)
        return set(
            template,
            lenses,
            (0, 1, atomic_positions, CellParametersCard("alat", cell_parameters.data)),
        )
    end
end # function update_structure

# This is a helper function and should not be exported.
function _preset(template::PWInput)
    @set! template.control = batchset(VerbositySetter(:high), template.control)
    @set! template.control.calculation = "scf"
    return batchset(CellParametersSetter(), template)
end # function _preset
function _preset(template::PhInput)
    lenses = @batchlens(begin
        _.inputph.verbosity  # Get the `template`'s `phonon.calculation` value
    end)
    # Set the `template`'s values with...
    template = set(template, lenses, ("high",))
end # function _preset

"""
    relay(from::PWInput, to::PhInput)

Relay shared information from a `PWInput` to a `PhInput`.

A `PWInput` before a `PhInput` has the information of `outdir` and `prefix`. They must keep the same in a
phonon calculation.
"""
function relay(from::PWInput, to::PhInput)
    lenses = @batchlens(begin
        _.outdir
        _.prefix
    end)
    return set(to, (@lens _.inputph) ∘ lenses, get(from, (@lens _.control) ∘ lenses))
end # function relay
"""
    relay(from::PhInput, to::Q2rInput)

Relay shared information from a `PhInput` to a `Q2rInput`.

A `PhInput` before a `Q2rInput` has the information of `fildyn`. It must keep the same in a q2r calculation.
"""
function relay(from::PhInput, to::Q2rInput)
    fildyn = @lens _.fildyn
    return set(to, (@lens _.input) ∘ fildyn, get(from, (@lens _.inputph) ∘ fildyn))
end # function relay
"""
    relay(from::Q2rInput, to::MatdynInput)

Relay shared information from a `Q2rInput` to a `MatdynInput`.

A `Q2rInput` before a `MatdynInput` has the information of `fildyn`, `flfrc` and `loto_2d`. They must keep the same
in a matdyn calculation.
"""
function relay(from::Q2rInput, to::MatdynInput)
    lenses = @batchlens(begin
        _.input.flfrc
        _.input.loto_2d
    end)
    return set(to, lenses, get(from, lenses))
end # function relay
"""
    relay(from::PhInput, to::DynmatInput)

Relay shared information from a `PhInput` to a `DynmatInput`.

A `PhInput` before a `DynmatInput` has the information of `asr`, `fildyn` and `amass`. They must keep the same
in a dynmat calculation.
"""
function relay(from::PhInput, to::DynmatInput)
    lenses = @batchlens(begin
        #_.asr   #TODO
        _.fildyn
        _.amass
    end)
    return set(to, (@lens _.input) ∘ lenses, get(from, (@lens _.inputph) ∘ lenses))
end # function relay

"""
    preprocess(step::Step{1}, inputs, outputs, template[, verbose::Bool = false])

Prepare input files of the first step of a phonon calculation.

# Arguments
- `step::Step{1}`: Denotes the first step of a phonon calculation.
- `inputs::AbstractVector{<:AbstractString}`: The input files of Quantum ESPRESSO of a phonon calculation. If not exist,
   they will be automatically created.
- `outputs::AbstractVector{<:AbstractString}`: The output files of Quantum ESPRESSO of a vc-relax calculation.
- `template::PWInput`:
- `verbose::Bool = false`: control the format of input files, verbose or not.
"""
function preprocess(
    ::Step{1},
    inputs::AbstractVector{<:AbstractString},
    outputs::AbstractVector{<:AbstractString},
    template::PWInput,
    verbose::Bool = false,
)
    # Check parameters
    @assert(
        length(inputs) == length(outputs),
        "The inputs, outputs must be the same length!"
    )
    template = _preset(template)
    for (input, output) in zip(inputs, outputs)
        # Get a new `object` from the `template`, with its `alat` and `pressure` changed
        object = update_structure(output, template)
        write(input, to_qe(object, verbose = verbose))  # Write the `object` to a Quantum ESPRESSO input file
    end
    return
end # function preprocess
function preprocess(
    ::Step{2},
    phonon_inputs::AbstractVector{<:AbstractString},
    pwscf_inputs::AbstractVector{<:AbstractString},
    template::PhInput,
    verbose::Bool = false,
)
    # Check parameters
    @assert(
        length(phonon_inputs) == length(pwscf_inputs),
        "The PWscf and the PHonon inputs files must have the same length!",
    )
    template = _preset(template)
    for (phonon_input, pwscf_input) in zip(phonon_inputs, pwscf_inputs)
        object = open(pwscf_input, "r") do io
            parse(PWInput, read(io, String))
        end
        template = relay(object, template)
        write(phonon_input, to_qe(template, verbose = verbose))
    end
    return
end # function preprocess
function preprocess(
    ::Step{3},
    q2r_inputs::AbstractVector{<:AbstractString},
    phonon_inputs::AbstractVector{<:AbstractString},
    template::Q2rInput,
    verbose::Bool = false,
)
    # Check parameters
    @assert(
        length(q2r_inputs) == length(phonon_inputs),
        "The phonon and the q2r inputs files must have the same length!",
    )
    for (q2r_input, phonon_input) in zip(q2r_inputs, phonon_inputs)
        object = open(phonon_input, "r") do io
            parse(PhInput, read(io, String))
        end
        template = relay(object, template)
        write(q2r_input, to_qe(template, verbose = verbose))
    end
    return
end # function preprocess
function preprocess(
    ::Step{4},
    matdyn_inputs::AbstractVector{<:AbstractString},
    q2r_inputs::AbstractVector{<:AbstractString},
    template::MatdynInput,
    verbose::Bool = false,
)
    # Check parameters
    @assert(
        length(matdyn_inputs) == length(q2r_inputs),
        "The q2r and the matdyn inputs files must have the same length!",
    )
    for (matdyn_input, q2r_input) in zip(matdyn_inputs, q2r_inputs)
        object = open(q2r_input, "r") do io
            parse(Q2rInput, read(io, String))
        end
        template = relay(object, template)
        if isfile(template.input.flfrq)
            @set! template.input.flfrq *= if template.input.dos == true  # Phonon DOS calculation
                "_dos"  # Append extension `"_dos` to `template.input.flfrq`
            else  # Phonon dispersion-relation calculation
                "_disp"  # Append extension `"_disp` to `template.input.flfrq`
            end
        end
        write(matdyn_input, to_qe(template, verbose = verbose))
    end
    return
end # function preprocess
function preprocess(
    ::Step{5},
    dynmat_inputs::AbstractVector{<:AbstractString},
    phonon_inputs::AbstractVector{<:AbstractString},
    template::DynmatInput,
    verbose::Bool = false,
)
    # Check parameters
    @assert(
        length(dynmat_inputs) == length(phonon_inputs),
        "The dynmat and the phonon inputs files must have the same length!",
    )
    for (dynmat_input, phonon_input) in zip(dynmat_inputs, phonon_inputs)
        object = open(phonon_input, "r") do io
            parse(PhInput, read(io, String))
        end
        template = relay(object, template)
        write(dynmat_input, to_qe(template, verbose = verbose))
    end
    return
end # function preprocess

#???


end
