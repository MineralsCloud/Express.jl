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

function set_press_vol(
    template::Input,
    pressure,
    eos::EquationOfState;
    minscale = eps(),
    maxscale = 1.3,
)
    @assert minscale > zero(minscale)  # No negative volume
)::Input
    volume = findvolume(eos(Pressure()), pressure, (minscale, maxscale) .* eos.v0)
    return _set_press_vol(template, pressure, volume)
end # function set_press_vol

const ALLOWED_CALCULATIONS = Union{SelfConsistentField,VariableCellOptimization}

function (step::Step{<:ALLOWED_CALCULATIONS,Prepare{:input}})(
    template::Input,
    pressure::Number,
    trial_eos::EquationOfState;
    minscale = eps(),
    maxscale = 1.3,
)
    template = preset(step, template)
    return set_press_vol(
        template,
        pressure,
        trial_eos;
        minscale = minscale,
        maxscale = maxscale,
    )
end
function (step::Step{<:ALLOWED_CALCULATIONS,Prepare{:input}})(
    input,
    template::Input,
    pressure::Number,
    trial_eos::EquationOfState;
    dry_run = false,
    kwargs...,
)
    object = step(template, pressure, trial_eos; kwargs...)
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
    return
end
function (step::Step{<:ALLOWED_CALCULATIONS,Prepare{:input}})(
    templates,
    pressures,
    trial_eos::EquationOfState;
    kwargs...,
)
    alert_pressures(pressures)
    return map(templates, pressures) do template, pressure  # `map` will check size mismatch
        step(template, pressure, trial_eos; kwargs...)
    end
end
function (step::Step{<:ALLOWED_CALCULATIONS,Prepare{:input}})(
    inputs,
    templates,
    pressures,
    trial_eos::EquationOfState,
    args...;
    dry_run = false,
    kwargs...,
)
    alert_pressures(pressures)
    map(inputs, templates, pressures) do input, template, pressure  # `map` will check size mismatch
        step(input, template, pressures, trial_eos, args...; dry_run = dry_run, kwargs...)
    end
    return
end
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

function (::Step{<:ALLOWED_CALCULATIONS,Launch{:job}})(
    outputs,
    inputs,
    n,
    bin;
    dry_run = false,
)
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
function (step::Step{T,Launch{:job}})(path::AbstractString) where {T<:ALLOWED_CALCULATIONS}
    settings = load_settings(path)
    inputs =
        @. settings.dirs * '/' * (T <: SelfConsistentField ? "scf" : "vc-relax") * ".in"
    outputs = map(Base.Fix2(replace, ".in" => ".out"), inputs)
    return step(outputs, inputs, settings.manager.np, settings.bin)
end

function (step::Step{<:ALLOWED_CALCULATIONS,Analyse{:output}})(
    outputs,
    trial_eos::EquationOfState,
    fit_e::Bool = true,
)
    results = map(outputs) do output
        analyse(step, read(output, String))  # volume => energy
    end
    if length(results) <= 5
        @info "pressures <= 5 may give unreliable results, run more if possible!"
    end
    if fit_e
        return lsqfit(trial_eos(Energy()), first.(results), last.(results))
    else
        return lsqfit(trial_eos(Pressure()), first.(results), last.(results))
    end
end
function (step::Step{VariableCellOptimization,Analyse{:output}})(output, template::Input)
    cell = open(output, "r") do io
        str = read(io, String)
        parsecell(str)
    end
    return set_structure(template, cell...)
end
function (step::Step{VariableCellOptimization,Analyse{:output}})(outputs, templates)
    return map(templates, outputs) do template, output  # `map` will check size mismatch
        step(output, template)
    end
end
function (step::Step{SelfConsistentField,Analyse{:output}})(path)
    settings = load_settings(path)
    inputs = settings.dirs .* "/scf.in"
    outputs = map(Base.Fix2(replace, ".in" => ".out"), inputs)
    return step(outputs, settings.trial_eos)
end
function (step::Step{VariableCellOptimization,Analyse{:output}})(path)
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
    if length(pressures) <= 5
        @info "pressures <= 5 may give unreliable results, consider more if possible!"
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

function set_structure end

function parsecell end

include("QuantumESPRESSO.jl")

end
