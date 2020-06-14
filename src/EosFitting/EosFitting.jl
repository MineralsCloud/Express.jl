"""
# module EosFitting



# Examples

```jldoctest
julia>
```
"""
module EosFitting

using AbInitioSoftwareBase.Inputs: Input, inputstring
using EquationsOfState.Collections: Pressure, Energy, EquationOfState
using EquationsOfState.NonlinearFitting: lsqfit
using EquationsOfState.Find: findvolume
using QuantumESPRESSO.CLI: pwcmd

using ..Express:
    Step,
    SelfConsistentField,
    VariableCellOptimization,
    Prepare,
    Launch,
    Analyse,
    PREPARE_POTENTIAL,
    PREPARE_INPUT,
    LAUNCH_JOB,
    ANALYSE_OUTPUT,
    load_settings,
    calculationtype,
    actiontype
using ..Jobs: nprocs_task, launchjob
using ..CLI: mpicmd

import ..Express

export Step,
    SelfConsistentField,
    VariableCellOptimization,
    PREPARE_POTENTIAL,
    PREPARE_INPUT,
    LAUNCH_JOB,
    ANALYSE_OUTPUT,
    load_settings,
    set_press_vol,
    inputstring,
    calculationtype,
    actiontype

mutable struct WorkingFiles
    pending
    running
    finished
end

function set_press_vol(
    template,
    pressure,
    eos::EquationOfState;
    minscale = eps(),
    maxscale = 1.3,
)
    @assert minscale > zero(minscale)  # No negative volume
    volume = findvolume(eos(Pressure()), pressure, (minscale, maxscale) .* eos.v0)
    return _set_press_vol(template, pressure, volume)
end # function set_press_vol

const ALLOWED_CALCULATIONS = Union{SelfConsistentField,VariableCellOptimization}

function (step::Step{<:ALLOWED_CALCULATIONS,Prepare{:input}})(
    f::Function,
    template::Input,
    pressure::Number,
    trial_eos::EquationOfState,
    args...;
    minscale = eps(),
    maxscale = 1.3,
    kwargs...,
)
    template = f(step, template, pressure, trial_eos, args...; kwargs...)  # A callback
    return set_press_vol(
        template,
        pressure,
        trial_eos;
        minscale = minscale,
        maxscale = maxscale,
    )
end
(step::Step{<:ALLOWED_CALCULATIONS,Prepare{:input}})(
    template::Input,
    pressure::Number,
    trial_eos::EquationOfState,
    args...;
    kwargs...,
) = step(preset, template, pressure, trial_eos, args...; kwargs...)
function (step::Step{<:ALLOWED_CALCULATIONS,Prepare{:input}})(
    f::Function,
    templates,
    pressures,
    trial_eos::EquationOfState,
    args...;
    kwargs...,
)
    return map(templates, pressures) do template, pressure  # `map` will check size mismatch
        step(f, template, pressure, trial_eos, args...; kwargs...)
    end
end
(step::Step{<:ALLOWED_CALCULATIONS,Prepare{:input}})(
    templates,
    pressures,
    trial_eos::EquationOfState,
    args...;
    kwargs...,
) = step(preset, templates, pressures, trial_eos, args...; kwargs...)
(step::Step{<:ALLOWED_CALCULATIONS,Prepare{:input}})(
    f::Function,
    template::Input,
    pressures,
    trial_eos::EquationOfState,
    args...;
    kwargs...,
) = step(f, fill(template, size(pressures)), pressures, trial_eos, args...; kwargs...)
(step::Step{<:ALLOWED_CALCULATIONS,Prepare{:input}})(
    template::Input,
    pressures,
    trial_eos::EquationOfState,
    args...;
    kwargs...,
) = step(fill(template, size(pressures)), pressures, trial_eos, args...; kwargs...)
function (step::Step{<:ALLOWED_CALCULATIONS,Prepare{:input}})(
    f::Function,
    inputs,
    templates,
    pressures,
    trial_eos::EquationOfState,
    args...;
    dry_run = false,
    kwargs...,
)
    objects = step(f, templates, pressures, trial_eos, args...; kwargs...)
    map(inputs, objects) do input, object  # `map` will check size mismatch
        if dry_run
            if isfile(input)
                @warn "file `$input` will be overwritten!"
            else
                @warn "file `$input` will be created!"
            end
            print(inputstring(object))
        else
            mkpath(dirname(input))
            open(input, "w") do io
                write(io, inputstring(object))
            end
        end
    end
    return
end
(step::Step{<:ALLOWED_CALCULATIONS,Prepare{:input}})(
    inputs,
    templates,
    pressures,
    trial_eos::EquationOfState,
    args...;
    dry_run = false,
    kwargs...,
) = step(preset, inputs, templates, pressures, trial_eos, args...; kwargs...)
function (step::Step{SelfConsistentField,Prepare{:input}})(path; kwargs...)
    settings = load_settings(path)
    inputs = settings.dirs .* "/scf.in"
    return step(
        inputs,
        settings.template,
        settings.pressures,
        settings.trial_eos;
        kwargs...,
    )
end
function (step::Step{VariableCellOptimization,Prepare{:input}})(path; kwargs...)
    settings = load_settings(path)
    inputs = settings.dirs .* "/vc-relax.in"
    new_eos = SelfConsistentField()(ANALYSE_OUTPUT)(path)
    return step(inputs, settings.template, settings.pressures, new_eos; kwargs...)
end

function (::Step{T,Launch{:job}})(outputs, inputs, n, bin; dry_run = false) where {T}
    # `map` guarantees they are of the same size, no need to check.
    n = nprocs_task(n, length(inputs))
    cmds = map(inputs, outputs) do input, output  # A vector of `Cmd`s
        _generate_cmds(n, input, output, bin)
    end
    if dry_run
        return cmds
    else
        return launchjob(cmds)
    end
end
function (step::Step{T,Launch{:job}})(path::AbstractString) where {T}
    settings = load_settings(path)
    inputs =
        @. settings.dirs * '/' * (T <: SelfConsistentField ? "scf" : "vc-relax") * ".in"
    outputs = map(Base.Fix2(replace, ".in" => ".out"), inputs)
    return step(outputs, inputs, settings.manager.np, settings.bin)
end

function (step::Step{T,Analyse{:output}})(outputs, trial_eos) where {T}
    results = map(outputs) do output
        str = read(output, String)
        analyse(step, str)  # volume => energy
    end
    return lsqfit(trial_eos(Energy()), first.(results), last.(results))
end # function postprocess
function (step::Step{SelfConsistentField,Analyse{:output}})(path::AbstractString)
    settings = load_settings(path)
    inputs = settings.dirs .* "/scf.in"
    outputs = map(Base.Fix2(replace, ".in" => ".out"), inputs)
    return step(outputs, settings.trial_eos)
end
function (step::Step{VariableCellOptimization,Analyse{:output}})(path::AbstractString)
    settings = load_settings(path)
    inputs = settings.dirs .* "/vc-relax.in"
    outputs = map(Base.Fix2(replace, ".in" => ".out"), inputs)
    new_eos = SelfConsistentField()(ANALYSE_OUTPUT)(path)
    return step(outputs, new_eos)
end # function preprocess

function (step::Step{SelfConsistentField,Prepare{:potential}})(template)
    required = getpotentials(template)
    path = getpotentialdir(template)
    return map(required) do potential
        download_potential(potential, path)
    end
end

# function (::T)(
#     outputs,
#     inputs,
#     template,
#     pressures,
#     trial_eos,
#     environment,
#     cmd,
# ) where {T<:Union{SelfConsistentField,VariableCellOptimization}}
#     Step{typeof(T),Prepare{:input}}(inputs, template, pressures, trial_eos)
#     Step{typeof(T),Launch{:job}}(outputs, inputs, environment, cmd)
#     Step{typeof(T),Analyse}(outputs, trial_eos)
# end

function Express._check_settings(settings)
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

# _generate_cmds(n, input, output, env::DockerEnvironment) = join(
#     [
#         "sh -c 'mpiexec --mca btl_vader_single_copy_mechanism none -np $n",
#         string('"', pwcmd(bin = env.bin).exec..., '"'),
#         "-inp \"$input\"'",
#     ],
#     " ",
# )
_generate_cmds(n, input, output, bin) =
    pipeline(mpicmd(n, pwcmd(bin = bin)), stdin = input, stdout = output)

function alert_pressures(pressures)
    if length(pressures) <= 6
        @info "pressures <= 6 may give unreliable results, consider more if possible!"
    end
    if minimum(pressures) >= zero(eltype(pressures))
        @warn "for better fitting, we need at least 1 negative pressure!"
    end
end # function alert_pressures

function _check_software_settings end

function _set_press_vol end

function preset end

function getpotentials end

function getpotentialdir end

function download_potential end

function analyse end

include("QuantumESPRESSO.jl")

end
