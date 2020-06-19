"""
# module Phonon



# Examples

```jldoctest
julia>
```
"""
module Phonon

using ConstructionBase: setproperties
using Kaleido: @batchlens
using QuantumESPRESSO.Inputs: InputFile, qestring
using QuantumESPRESSO.Inputs.PWscf:
    AtomicPositionsCard, CellParametersCard, PWInput, optconvert
using QuantumESPRESSO.Inputs.PHonon: PhInput, Q2rInput, MatdynInput, DynmatInput
using QuantumESPRESSO.Outputs: OutputFile
using QuantumESPRESSO.Outputs.PWscf: parsefinal
using QuantumESPRESSO.Setters: CellParametersSetter, VerbositySetter
using Setfield: get, set, @lens, @set!

using ..Express: Calculation, SelfConsistentField, PREPARE_INPUT, LAUNCH_JOB, ANALYSE_OUTPUT

export set_structure, relay, preprocess

struct DensityFunctionalPertubation{T} <: Calculation{T} end
struct ReciprocalToReal{T} <: Calculation{T} end
struct PhononDensityOfStates{T} <: Calculation{T} end
struct PhononDispersion{T} <: Calculation{T} end

"""
    set_structure(output::AbstractString, template::PWInput)

Read structure information from `output`, and update related fields of `template`.
"""
function set_structure(output::AbstractString, template::PWInput)
    str = read(OutputFile(output))
    cell = parsefinal(CellParametersCard{Float64}, str)
    return setproperties(
        template,
        atomic_positions = parsefinal(AtomicPositionsCard, str),
        cell_parameters = CellParametersCard(cell.data / template.system.celldm[1], "alat"), # The result of `parsefinal` must be a `CellParametersCard` with `"bohr"` or `"angstrom"` option, convert it to "bohr" by default
    )
end # function set_structure

# This is a helper function and should not be exported.
function preset(template::PWInput)
    @set! template.control = set(template.control, VerbositySetter(:high))
    @set! template.control.calculation = "scf"
    template = set(template, CellParametersSetter())
    @set! template.system.ibrav = 0
    return template
end # function preset
function preset(template::PhInput)
    lenses = @batchlens(begin
        _.inputph.verbosity  # Get the `template`'s `phonon.calculation` value
    end)
    # Set the `template`'s values with...
    template = set(template, lenses, ("high",))
end # function preset

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
"""
function (::SelfConsistentField{PREPARE_INPUT})(scf_inputs, vc_outputs, template::PWInput)
    template = preset(template)
    map(scf_inputs, vc_outputs) do input, output
        # Get a new `object` from the `template`, with its `alat` and `pressure` changed
        object = set_structure(output, template)
        write(InputFile(input), object)
        object
    end
end # function preprocess
function preprocess(
    ::DensityFunctionalPertubation{PREPARE_INPUT},
    phonon_inputs,
    pwscf_inputs,
    template::PhInput,
)
    template = preset(template)
    map(phonon_inputs, pwscf_inputs) do phonon_input, pwscf_input
        object = parse(PWInput, read(InputFile(pwscf_input)))
        write(InputFile(phonon_input), relay(object, template))
    end
    return
end # function preprocess
function preprocess(
    ::ReciprocalToReal{PREPARE_INPUT},
    q2r_inputs,
    phonon_inputs,
    template::Q2rInput,
)
    map(q2r_inputs, phonon_inputs) do q2r_input, phonon_input
        object = parse(PhInput, read(InputFile(phonon_input)))
        write(InputFile(q2r_input), relay(object, template))
    end
    return
end # function preprocess
function preprocess(
    ::PhononDispersion{PREPARE_INPUT},
    matdyn_inputs,
    q2r_inputs,
    template::MatdynInput,
)
    map(matdyn_inputs, q2r_inputs) do matdyn_input, q2r_input
        object = parse(Q2rInput, read(InputFile(q2r_input)))
        template = relay(object, template)
        if isfile(template.input.flfrq)
            @set! template.input.flfrq *= if template.input.dos == true  # Phonon DOS calculation
                "_dos"  # Append extension `"_dos` to `template.input.flfrq`
            else  # Phonon dispersion-relation calculation
                "_disp"  # Append extension `"_disp` to `template.input.flfrq`
            end
        end
        write(InputFile(matdyn_input), template)
    end
    return
end # function preprocess
function preprocess(
    ::PhononDispersion{PREPARE_INPUT},
    dynmat_inputs,
    phonon_inputs,
    template::DynmatInput,
)
    map(dynmat_inputs, phonon_inputs) do dynmat_input, phonon_input
        object = parse(PhInput, read(InputFile(phonon_input)))
        write(InputFile(dynmat_input), relay(object, template))
    end
    return
end # function preprocess

#???

end
