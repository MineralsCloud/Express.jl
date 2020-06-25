"""
# module EosFitting



# Examples

```jldoctest
julia>
```
"""
module EosFitting

using AbInitioSoftwareBase.Inputs: Input, inputstring
using AbInitioSoftwareBase.CLI: MpiCmd
using EquationsOfState.Collections: Pressure, Energy, EquationOfState
using EquationsOfState.NonlinearFitting: lsqfit
using EquationsOfState.Find: findvolume
using OptionalArgChecks: @argcheck

using ..Express: SelfConsistentField, VariableCellOptimization, load_settings
using ..Jobs: nprocs_task, launchjob

import ..Express

export SelfConsistentField,
    VariableCellOptimization,
    load_settings,
    set_press_vol,
    inputstring,
    preprocess,
    process,
    postprocess

const ALLOWED_CALCULATIONS = Union{SelfConsistentField,VariableCellOptimization}

function set_press_vol(
    template::Input,
    pressure,
    eos::EquationOfState;
    volume_scale = (eps(), 1.3),
)::Input
    ⋁, ⋀ = minimum(volume_scale), maximum(volume_scale)
    @argcheck ⋁ > zero(eltype(volume_scale))  # No negative volume
    volume = findvolume(eos(Pressure()), pressure, (⋁, ⋀) .* eos.v0)
    return _set_press_vol(template, pressure, volume)
end # function set_press_vol

function prep_input(
    calc::ALLOWED_CALCULATIONS,
    template::Input,
    pressure::Number,
    trial_eos::EquationOfState;
    volume_scale = (eps(), 1.3),
)
    return set_press_vol(_prep_input(calc, template), pressure, trial_eos, volume_scale)
end
function prep_input(
    calc::ALLOWED_CALCULATIONS,
    templates,
    pressures,
    trial_eos::EquationOfState;
    kwargs...,
)
    alert_pressures(pressures)
    return map(templates, pressures) do template, pressure  # `map` will check size mismatch
        prep_input(calc, template, pressure, trial_eos; kwargs...)
    end
end

function write_input(file, object::Input; dry_run = false)
    if dry_run
        if isfile(file)
            @warn "file `$file` will be overwritten!"
        else
            @warn "file `$file` will be created!"
        end
        print(inputstring(object))
    else
        mkpath(dirname(file))
        open(file, "w") do io
            write(io, inputstring(object))
        end
    end
    return
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
        object = prep_input(calc, template, pressures, trial_eos; kwargs...)
        write_input(file, object; dry_run = dry_run)
    end
    return
end
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

function process(outputs, inputs, n, softwarecmd; dry_run = false, kwargs...)
    # `map` guarantees they are of the same size, no need to check.
    n = nprocs_task(n, length(inputs))
    cmds = map(inputs, outputs) do input, output  # A vector of `Cmd`s
        f = MpiCmd(n; kwargs...) ∘ softwarecmd
        f(stdin = input, stdout = output)
    end
    if dry_run
        @warn "the following commands will be run:"
        return cmds
    else
        return launchjob(cmds)
    end
end
function process(calc::T, path::AbstractString) where {T<:ALLOWED_CALCULATIONS}
    settings = load_settings(path)
    inputs =
        @. settings.dirs * '/' * (T <: SelfConsistentField ? "scf" : "vc-relax") * ".in"
    outputs = map(Base.Fix2(replace, ".in" => ".out"), inputs)
    return process(calc, outputs, inputs, settings.manager.np, settings.bin)
end

function postprocess(
    outputs,
    trial_eos::EquationOfState,
    fit_e::Bool = true,
)::EquationOfState
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
        @argcheck haskey(settings, key)
    end
    _check_software_settings(settings["qe"])
    @argcheck isdir(settings["dir"])
    @argcheck isfile(settings["template"])
    alert_pressures(settings["pressures"])
    map(("type", "parameters", "units")) do key
        @argcheck haskey(settings["trial_eos"], key)
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

function _prep_input end

function getpotentials end

function getpotentialdir end

function download_potential end

function analyse end

function set_structure end

function parsecell end

include("QuantumESPRESSO.jl")

end
