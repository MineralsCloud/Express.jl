"""
# module Phonon



# Examples

```jldoctest
julia>
```
"""
module Phonon

using AbInitioSoftwareBase: loadfile
using AbInitioSoftwareBase.Inputs: Input, inputstring, writeinput

using ..Express:
    SelfConsistentField,
    DfptMethod,
    ForceConstant
using ..EosFitting: _check_software_settings
import ..Express

export DfptMethod, ForceConstant, prep_input, prepare, process, load_settings

function relay end

function prep_input end

function (step::Step{VariableCellOptimization,Action{:set_structure}})(
    outputs,
    template::Input,
)
    map(outputs) do output
        cell = open(output, "r") do io
            str = read(io, String)
            parsecell(str)
        end
        if any(x === nothing for x in cell)
            return
        else
            return _set_structure(template, cell...)
        end
    end
end
function (step::Step{VariableCellOptimization,Action{:set_structure}})(configfile)
    settings = load_settings(configfile)
    inputs = settings.dirs .* "/vc-relax.in"
    outputs = map(Base.Fix2(replace, ".in" => ".out"), inputs)
    new_inputs = settings.dirs .* "/new.in"
    for (st, input) in zip(step(outputs, settings.template), new_inputs)
        if st !== nothing
            write(input, inputstring(st))
        end
    end
end

function prepare(::DfptMethod, inputs, template::Input, args...; dry_run = false)
    map(inputs) do input
        writeinput(input, prep_input(DfptMethod(), template, args...), dry_run)
    end
    return
end
function prepare(calc::DfptMethod, path; kwargs...)
    settings = load_settings(path)
    inputs = settings.dirs .* "/ph.in"
    return prepare(calc, inputs, settings.template[2], settings.template[1]; kwargs...)
end
function prepare(::ForceConstant, inputs, template::Input, args...; dry_run = false)
    map(inputs) do input
        Step(ForceConstant(), PREPARE_INPUT)(input, template, args...)
    end
    return
end
function prepare(calc::ForceConstant, path; kwargs...)
    settings = load_settings(path)
    inputs = settings.dirs .* "/q2r.in"
    return prepare(calc, inputs, settings.template[2], settings.template[1]; kwargs...)
end

function process(
    calc,
    outputs,
    inputs,
    n,
    softwarecmd;
    dry_run = false,
    kwargs...,
)
    if !dry_run
        Step(calc, LAUNCH_JOB)(outputs, inputs, n, softwarecmd)
    end
end
function process(calc::T, path::AbstractString; kwargs...) where {T}
    settings = load_settings(path)
    inputs =
        @. settings.dirs * '/' * (T <: SelfConsistentField ? "scf" : "ph") * ".in"
    outputs = map(Base.Fix2(replace, ".in" => ".out"), inputs)
    return process(calc, outputs, inputs, settings.manager.np, settings.bin; kwargs...)
end

# function (::Step{ForceConstant,Prepare{:input}})(
#     q2r_inputs,
#     phonon_inputs,
#     template::Q2rInput,
# )
#     map(q2r_inputs, phonon_inputs) do q2r_input, phonon_input
#         object = parse(PhInput, read(InputFile(phonon_input)))
#         write(InputFile(q2r_input), relay(object, template))
#     end
#     return
# end

# function (::Step{PhononDispersion,Prepare{:input}})(
#     matdyn_inputs,
#     q2r_inputs,
#     template::MatdynInput,
# )
#     map(matdyn_inputs, q2r_inputs) do matdyn_input, q2r_input
#         object = parse(Q2rInput, read(InputFile(q2r_input)))
#         template = relay(object, template)
#         if isfile(template.input.flfrq)
#             @set! template.input.flfrq *= if template.input.dos == true  # Phonon DOS calculation
#                 "_dos"  # Append extension `"_dos` to `template.input.flfrq`
#             else  # Phonon dispersion-relation calculation
#                 "_disp"  # Append extension `"_disp` to `template.input.flfrq`
#             end
#         end
#         write(InputFile(matdyn_input), template)
#     end
#     return
# end
# function (::Step{PhononDispersion,Prepare{:input}})(
#     dynmat_inputs,
#     phonon_inputs,
#     template::DynmatInput,
# )
#     map(dynmat_inputs, phonon_inputs) do dynmat_input, phonon_input
#         object = parse(PhInput, read(InputFile(phonon_input)))
#         write(InputFile(dynmat_input), relay(object, template))
#     end
#     return
# end

function _expand_settings end

function _check_settings(settings)
    map(("template", "pressures", "dir")) do key
        @assert haskey(settings, key)
    end
    @assert isdir(settings["dir"])
    @assert all(isfile.(settings["template"]))
end # function _check_settings

function load_settings(path)
    settings = loadfile(path)
    _check_settings(settings)  # Errors will be thrown if exist
    return _expand_settings(settings)
end # function load_settings

include("QuantumESPRESSO.jl")

end
