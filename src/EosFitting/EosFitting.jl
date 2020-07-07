"""
# module EosFitting



# Examples

```jldoctest
julia>
```
"""
module EosFitting

using AbInitioSoftwareBase: load
using AbInitioSoftwareBase.Inputs: Input, inputstring, write_input
using AbInitioSoftwareBase.CLI: MpiCmd
using Dates: DateTime, now, format
using EquationsOfState.Collections: Pressure, Energy, EquationOfState
using EquationsOfState.NonlinearFitting: lsqfit
using EquationsOfState.Find: findvolume

using ..Express: SelfConsistentField, VariableCellOptimization, Calculation
using ..Jobs: div_nprocs, launchjob

import ..Express

export SelfConsistentField,
    VariableCellOptimization,
    load_settings,
    set_press_vol,
    inputstring,
    preprocess,
    process,
    postprocess,
    set_structure

const ALLOWED_CALCULATIONS = Union{SelfConsistentField,VariableCellOptimization}

@enum StepStatus::Bool begin
    SUCCEEDED
    FAILED
end

function set_press_vol(
    template::Input,
    pressure,
    eos::EquationOfState;
    volume_scale = (eps(), 1.3),
)::Input
    ⋁, ⋀ = minimum(volume_scale), maximum(volume_scale)
    @assert ⋁ > zero(eltype(volume_scale))  # No negative volume
    volume = findvolume(eos(Pressure()), pressure, (⋁, ⋀) .* eos.v0)
    return _set_press_vol(template, pressure, volume)
end # function set_press_vol

function prep_input(
    calc::ALLOWED_CALCULATIONS,
    template::Input,
    pressure::Number,
    trial_eos::EquationOfState;
    kwargs...,
)
    return set_press_vol(_prep_input(calc, template), pressure, trial_eos; kwargs...)
end

function preprocess(
    calc::ALLOWED_CALCULATIONS,
    files,
    templates,
    pressures,
    trial_eos::EquationOfState;
    dry_run = false,
    kwargs...,
)
    alert_pressures(pressures)
    map(files, templates, pressures) do file, template, pressure  # `map` will check size mismatch
        object = prep_input(calc, template, pressure, trial_eos; kwargs...)
        write_input(file, object, dry_run)
    end
    STEP_TRACKER[calc isa SelfConsistentField ? 1 : 4] =
        Context(files, nothing, SUCCEEDED, now(), calc)
    return
end
preprocess(
    calc::ALLOWED_CALCULATIONS,
    files,
    template::Input,
    pressures,
    trial_eos::EquationOfState;
    kwargs...,
) = preprocess(calc, files, fill(template, size(files)), pressures, trial_eos; kwargs...)
function preprocess(calc::SelfConsistentField, path; kwargs...)
    settings = load_settings(path)
    inputs = settings.dirs .* "/scf.in"
    return preprocess(
        calc,
        inputs,
        settings.template,
        settings.pressures,
        settings.trial_eos;
        kwargs...,
    )
end
function preprocess(calc::VariableCellOptimization, path; kwargs...)
    settings = load_settings(path)
    inputs = settings.dirs .* "/vc-relax.in"
    new_eos = preprocess(SelfConsistentField(), path)
    return preprocess(
        calc,
        inputs,
        settings.template,
        settings.pressures,
        new_eos;
        kwargs...,
    )
end

function process(
    calc::ALLOWED_CALCULATIONS,
    outputs,
    inputs,
    n,
    softwarecmd;
    dry_run = false,
    kwargs...,
)
    # `map` guarantees they are of the same size, no need to check.
    n = div_nprocs(n, length(inputs))
    cmds = map(inputs, outputs) do input, output  # A vector of `Cmd`s
        f = MpiCmd(n; kwargs...) ∘ softwarecmd
        f(stdin = input, stdout = output)
    end
    if dry_run
        @warn "the following commands will be run:"
        return cmds
    else
        STEP_TRACKER[calc isa SelfConsistentField ? 2 : 5] =
            Context(inputs, outputs, SUCCEEDED, now(), calc)
        return launchjob(cmds)
    end
end
function process(calc::T, path::AbstractString; kwargs...) where {T<:ALLOWED_CALCULATIONS}
    settings = load_settings(path)
    inputs =
        @. settings.dirs * '/' * (T <: SelfConsistentField ? "scf" : "vc-relax") * ".in"
    outputs = map(Base.Fix2(replace, ".in" => ".out"), inputs)
    return process(calc, outputs, inputs, settings.manager.np, settings.bin, kwargs...)
end

function postprocess(
    calc::ALLOWED_CALCULATIONS,
    outputs,
    trial_eos::EquationOfState,
    fit_e::Bool = true,
)::EquationOfState
    results = map(outputs) do output
        analyse(calc, read(output, String))  # volume => energy
    end
    if length(results) <= 5
        @info "pressures <= 5 may give unreliable results, run more if possible!"
    end
    STEP_TRACKER[calc isa SelfConsistentField ? 3 : 6] =
        Context(nothing, outputs, SUCCEEDED, now(), calc)
    if fit_e
        return lsqfit(trial_eos(Energy()), first.(results), last.(results))
    else
        return lsqfit(trial_eos(Pressure()), first.(results), last.(results))
    end
end
function postprocess(calc::SelfConsistentField, path)
    settings = load_settings(path)
    inputs = settings.dirs .* "/scf.in"
    outputs = map(Base.Fix2(replace, ".in" => ".out"), inputs)
    return postprocess(calc, outputs, settings.trial_eos)
end
function postprocess(::VariableCellOptimization, path)
    settings = load_settings(path)
    inputs = settings.dirs .* "/vc-relax.in"
    outputs = map(Base.Fix2(replace, ".in" => ".out"), inputs)
    new_eos = postprocess(SelfConsistentField(), path)
    return postprocess(VariableCellOptimization(), outputs, new_eos)
end

function set_structure(::VariableCellOptimization, output, template::Input)
    cell = open(output, "r") do io
        str = read(io, String)
        parsecell(str)
    end
    return set_structure(template, cell...)
end
function set_structure(::VariableCellOptimization, outputs, templates)
    return map(templates, outputs) do template, output  # `map` will check size mismatch
        step(output, template)
    end
end

function prep_potential(template)
    required = getpotentials(template)
    path = getpotentialdir(template)
    return map(required) do potential
        download_potential(potential, path)
    end
end

mutable struct Context
    inputs
    outputs
    status::StepStatus
    time::DateTime
    calc::Calculation
end

STEP_TRACKER = [
    Context(nothing, nothing, FAILED, now(), SelfConsistentField()),
    Context(nothing, nothing, FAILED, now(), SelfConsistentField()),
    Context(nothing, nothing, FAILED, now(), SelfConsistentField()),
    Context(nothing, nothing, FAILED, now(), VariableCellOptimization()),
    Context(nothing, nothing, FAILED, now(), VariableCellOptimization()),
    Context(nothing, nothing, FAILED, now(), VariableCellOptimization()),
]

function Base.show(io::IO, x::Context)
    print(io, _emoji(x.status), ' ')
    printstyled(io, x.calc; bold = true)
    printstyled(io, " @ ", format(x.time, "Y/mm/dd H:M:S"); color = :light_black)
end # function Base.show

_emoji(step) = step == SUCCEEDED ? '✅' : '❌'

# _generate_cmds(n, input, output, env::DockerEnvironment) = join(
#     [
#         "sh -c 'mpiexec --mca btl_vader_single_copy_mechanism none -np $n",
#         string('"', pwcmd(bin = env.bin).exec..., '"'),
#         "-inp \"$input\"'",
#     ],
#     " ",
# )

function alert_pressures(pressures)
    if length(pressures) <= 5
        @info "pressures <= 5 may give unreliable results, consider more if possible!"
    end
    if minimum(pressures) >= zero(eltype(pressures))
        @warn "for better fitting, we need at least 1 negative pressure!"
    end
end # function alert_pressures

function _set_press_vol end

function _prep_input end

function getpotentials end

function getpotentialdir end

function download_potential end

function analyse end

function set_structure end

function parsecell end

function _expand_settings end

function _check_software_settings end

function _check_settings(settings)
    map(("template", "pressures", "trial_eos", "dir")) do key
        @assert haskey(settings, key)
    end
    _check_software_settings(settings["qe"])
    @assert isdir(settings["dir"])
    @assert isfile(settings["template"])
    alert_pressures(settings["pressures"])
    map(("type", "parameters", "units")) do key
        @assert haskey(settings["trial_eos"], key)
    end
end # function _check_settings

function load_settings(path)
    settings = load(path)
    _check_settings(settings)  # Errors will be thrown if exist
    return _expand_settings(settings)
end # function load_settings

include("QuantumESPRESSO.jl")

end
