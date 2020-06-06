"""
# module EosFitting



# Examples

```jldoctest
julia>
```
"""
module EosFitting

using Crystallography: cellvolume
using EquationsOfState.Collections: Pressure, Energy, EquationOfState
using EquationsOfState.NonlinearFitting: lsqfit
using EquationsOfState.Find: findvolume
using Unitful: NoUnits
using UnitfulAtomic: bohr, Ry
using QuantumESPRESSO.CLI: pwcmd

using ..Express:
    Step,
    SelfConsistentField,
    VariableCellRelaxation,
    PrepareInput,
    LaunchJob,
    AnalyseOutput,
    load_settings,
    inputstring
using ..Jobs: nprocs_task, distribute_process
using ..Environments: DockerEnvironment, LocalEnvironment
using ..CLI: mpicmd

import ..Express

export Step,
    SelfConsistentField,
    VariableCellRelaxation,
    PrepareInput,
    LaunchJob,
    AnalyseOutput,
    load_settings,
    parse_template,
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

function (step::Step{T,PrepareInput})(
    template,
    pressures,
    trial_eos::EquationOfState;
    kwargs...,
) where {T}
    template = _set_boilerplate(T, template)
    return map(pressures) do pressure  # `map` will check size mismatch
        set_press_vol(template, pressure, trial_eos; kwargs...)  # Create a new `object` from `template`, with its `alat` and `pressure` changed
    end
end
function (step::Step{T,PrepareInput})(
    inputs,
    template,
    pressures,
    trial_eos::EquationOfState;
    dry_run = false,
    kwargs...,
) where {T}
    objects = step(template, pressures, trial_eos; kwargs...)
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
function (step::Step{T,PrepareInput})(path::AbstractString) where {T}
    settings = load_settings(path)
    return step(settings.inputs, settings.template, settings.pressures, settings.trial_eos)
end # function preprocess

function (::Step{T,LaunchJob})(outputs, inputs, environment; dry_run = false) where {T}
    # `map` guarantees they are of the same size, no need to check.
    n = nprocs_task(environment.n, length(inputs))
    cmds = map(inputs, outputs) do input, output  # A vector of `Cmd`s
        _generate_cmds(n, input, output, environment)
    end
    if dry_run
        return cmds
    else
        return distribute_process(cmds, environment)
    end
end
function (step::Step{T,LaunchJob})(path::AbstractString) where {T}
    settings = load_settings(path)
    outputs = map(Base.Fix2(replace, ".in" => ".out"), settings.inputs)
    return step(outputs, settings.inputs, settings.environment)
end

function (step::Step{T,AnalyseOutput})(outputs, trial_eos) where {T}
    results = map(outputs) do output
        s = read(output, String)
        parseenergies(step, s)
    end
    return lsqfit(trial_eos(Energy()), volumes(results), energies(results))
end # function postprocess
function (step::Step{T,AnalyseOutput})(path::AbstractString) where {T}
    settings = load_settings(path)
    outputs = map(Base.Fix2(replace, ".in" => ".out"), settings.inputs)
    return step(outputs, settings.trial_eos)
end # function (step::Step{T,AnalyseOutput})

# function (::T)(
#     outputs,
#     inputs,
#     template,
#     pressures,
#     trial_eos,
#     environment,
#     cmd,
# ) where {T<:Union{SelfConsistentField,VariableCellRelaxation}}
#     Step{typeof(T),PrepareInput}(inputs, template, pressures, trial_eos)
#     Step{typeof(T),LaunchJob}(outputs, inputs, environment, cmd)
#     Step{typeof(T),AnalyseOutput}(outputs, trial_eos)
# end

function Express._check_settings(settings)
    map(("template", "pressures", "trial_eos", "dir")) do key
        @assert haskey(settings, key)
    end
    _check_software_settings(settings["qe"])
    @assert isdir(settings["dir"])
    @assert isfile(settings["template"])
    if length(settings["pressures"]) <= 6
        @info "pressures less than 6 may give unreliable results, consider more if possible!"
    end
    map(("type", "parameters", "units")) do key
        @assert haskey(settings["trial_eos"], key)
    end
end # function _check_settings

_generate_cmds(n, input, output, env::DockerEnvironment) = join(
    [
        "sh -c 'mpiexec --mca btl_vader_single_copy_mechanism none -np $n",
        pwcmd(bin = env.bin).exec...,
        "-inp $input'",
    ],
    " ",
)
_generate_cmds(n, input, output, env::LocalEnvironment) =
    pipeline(mpicmd(n, pwcmd(bin = env.bin)), stdin = input, stdout = output)

mutable struct ContextManager
    environment::SimulationEnvironment
end

function parse_template end

function parseenergies end

function _set_boilerplate end

function _check_software_settings end

function _set_press_vol end

function volumes end

function energies end

include("QuantumESPRESSO.jl")

end
