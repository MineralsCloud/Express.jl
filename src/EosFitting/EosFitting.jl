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
    inputstring
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
    inputstring

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

function (step::Step{T,Prepare{:input}})(
    template::Input,
    pressures,
    trial_eos::EquationOfState;
    kwargs...,
) where {T}
    template = step(template)  # To be extended
    return step(fill(template, size(pressures)), pressures, trial_eos; kwargs...)
end
function (step::Step{T,Prepare{:input}})(
    templates,
    pressures,
    trial_eos::EquationOfState;
    kwargs...,
) where {T}
    alert_pressures(pressures)
    return map(templates, pressures) do template, pressure  # `map` will check size mismatch
        set_press_vol(template, pressure, trial_eos; kwargs...)  # Create a new `object` from `template`, with its `alat` and `pressure` changed
    end
end
function (step::Step{T,Prepare{:input}})(
    inputs,
    templates,
    pressures,
    trial_eos::EquationOfState;
    dry_run = false,
    kwargs...,
) where {T}
    objects = step(templates, pressures, trial_eos; kwargs...)
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
function (step::Step{SelfConsistentField,Prepare{:input}})(path::AbstractString)
    settings = load_settings(path)
    inputs = settings.dirs .* "/scf.in"
    return step(inputs, settings.template, settings.pressures, settings.trial_eos)
end # function preprocess
function (step::Step{VariableCellOptimization,Prepare{:input}})(path::AbstractString)
    settings = load_settings(path)
    inputs = settings.dirs .* "/vc-relax.in"
    new_eos = SelfConsistentField(ANALYSE_OUTPUT)(path, settings.trial_eos)
    return step(inputs, settings.template, settings.pressures, new_eos)
end # function preprocess

function (::Step{T,Launch{:job}})(outputs, inputs, environment; dry_run = false) where {T}
    # `map` guarantees they are of the same size, no need to check.
    n = nprocs_task(environment.n, length(inputs))
    cmds = map(inputs, outputs) do input, output  # A vector of `Cmd`s
        _generate_cmds(n, input, output, environment)
    end
    if dry_run
        return cmds
    else
        return launchjob(cmds, environment)
    end
end
function (step::Step{T,Launch{:job}})(path::AbstractString) where {T}
    settings = load_settings(path)
    inputs = @. settings.dirs * '/' * (T <: SelfConsistentField ? "scf" : "vc-relax") * ".in"
    outputs = map(Base.Fix2(replace, ".in" => ".out"), inputs)
    return step(outputs, inputs, settings.environment)
end

function (step::Step{T,Analyse{:output}})(outputs, trial_eos) where {T}
    strs = (read(output, String) for output in outputs)
    results = (step(str) for str in strs)  # [volume => energy]
    return lsqfit(trial_eos(Energy()), keys(results), values(results))
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
    new_eos = SelfConsistentField(ANALYSE_OUTPUT)(path, settings.trial_eos)
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
_generate_cmds(n, input, output, env) =
    pipeline(mpicmd(n, pwcmd(bin = env.bin)), stdin = input, stdout = output)

function alert_pressures(pressures)
    if length(pressures) <= 6
        @info "pressures <= 6 may give unreliable results, consider more if possible!"
    end
    if minimum(pressures) >= zero(eltype(pressures))
        @warn "for better fitting, we need at least 1 negative pressure!"
    end
end # function alert_pressures

mutable struct ContextManager
    environment::CalculationEnvironment
end

function _check_software_settings end

function _set_press_vol end

function getpotentials end

function getpotentialdir end

function download_potential end

include("QuantumESPRESSO.jl")

end
