module EosFitting

using AbInitioSoftwareBase: FilePath, load
using AbInitioSoftwareBase.Inputs: Input, inputstring, write_input
using AbInitioSoftwareBase.CLI: MpiExec
using Dates: now
using EquationsOfState.Collections: Pressure, Energy, EquationOfState
using EquationsOfState.NonlinearFitting: lsqfit
using EquationsOfState.Find: findvolume

using ..Express:
    Optimization,
    StructuralOptimization,
    VariableCellOptimization,
    SelfConsistentField,
    Calculation,
    Action,
    Step,
    Succeeded,
    Pending,
    Context,
    _emoji
using ..Jobs: div_nprocs, launchjob

import ..Express

export SelfConsistentField,
    StructuralOptimization,
    VariableCellOptimization,
    load_settings,
    set_pressure_volume,
    inputstring,
    prepare,
    process,
    postprocess,
    set_structure

const ScfOrOptim = Union{SelfConsistentField,Optimization}
const PREPARE_INPUT = Action{:prepare_input}()
const LAUNCH_JOB = Action{:launch_job}()
const WRITE_INPUT = Action{:write_input}()
const FIT_EOS = Action{:fit_eos}()
const SET_STRUCTURE = Action{:set_structure}()

"""
    set_pressure_volume(template::Input, pressure, eos::EquationOfState; volume_scale = (eps(), 1.3))

Set the volume of `template` at a `pressure` according to `eos`.

The `volume_scale` gives a trial of the minimum and maximum scales for the `eos`. It
times the zero-pressure volume of the `eos` will be the trial volumes.
"""
function set_pressure_volume(
    template::Input,
    pressure,
    eos::EquationOfState;
    volume_scale = (eps(), 1.3),
)::Input
    @assert minimum(volume_scale) > zero(eltype(volume_scale))  # No negative volume
    volume = findvolume(eos(Pressure()), pressure, extrema(volume_scale) .* eos.v0)
    return _set_pressure_volume(template, pressure, volume)
end # function set_pressure_volume

"""
    (step::Step{<:Union{SelfConsistentField,Optimization},Action{:prepare_input}})(template::Input, pressure, trial_eos::EquationOfState; kwargs...,)

Generate input files from a given `template`, `pressure` and `trial_eos`, with some preset values.

See also: [`set_pressure_volume`](@ref)
"""
function (step::Step{<:ScfOrOptim,Action{:prepare_input}})(
    template::Input,
    pressure::Number,
    trial_eos::EquationOfState;
    kwargs...,
)
    return set_pressure_volume(step(template), pressure, trial_eos; kwargs...)
end

function (step::Step{T,Action{:write_input}})(
    file,
    template::Input,
    pressure::Number,
    trial_eos::EquationOfState;
    dry_run = false,
    kwargs...,
) where {T<:ScfOrOptim}
    object = Step(T(), PREPARE_INPUT)(template, pressure, trial_eos; kwargs...)
    write_input(file, object, dry_run)
    return  # Always return `nothing` no matter `dry_run` is `true` or `false`.
end

"""
    preprocess(calc::Union{SelfConsistentField,VariableCellOptimization}, files, template::Input, pressures, trial_eos::EquationOfState; kwargs...)
    preprocess(calc::Union{SelfConsistentField,VariableCellOptimization}, files, templates, pressures, trial_eos::EquationOfState; dry_run, kwargs...)

Prepare the input `files` from a certain `template` / a series of `templates` at `pressures` from a `trial_eos`.

Set `dry_run = true` to see what will happen instead of actual happening.
"""
function prepare(
    calc::ScfOrOptim,
    files,
    templates,
    pressures,
    trial_eos::EquationOfState;
    dry_run = false,
    kwargs...,
)
    alert_pressures(pressures)
    step = Step(calc, WRITE_INPUT)
    for (file, template, pressure) in zip(files, templates, pressures)
        step(file, template, pressure, trial_eos; kwargs...)
    end
    STEP_TRACKER[calc isa SelfConsistentField ? 1 : 4] =
        Context(files, nothing, Succeeded(), now(), Step(calc, PREPARE_INPUT))
    return
end
prepare(
    calc::ScfOrOptim,
    files,
    template::Input,
    pressures,
    trial_eos::EquationOfState;
    kwargs...,
) = prepare(calc, files, fill(template, size(files)), pressures, trial_eos; kwargs...)
"""
    preprocess(calc::Union{SelfConsistentField,VariableCellOptimization}, configfile; kwargs...)

Do the same thing of `preprocess`, but from a configuration file.
"""
function prepare(calc::SelfConsistentField, configfile; kwargs...)
    settings = load_settings(configfile)
    inputs = settings.dirs .* "/scf.in"
    return prepare(
        calc,
        inputs,
        settings.template,
        settings.pressures,
        settings.trial_eos;
        kwargs...,
    )
end
function prepare(calc::VariableCellOptimization, configfile; kwargs...)
    settings = load_settings(configfile)
    inputs = settings.dirs .* "/vc-relax.in"
    new_eos = postprocess(SelfConsistentField(), configfile)
    return prepare(
        calc,
        inputs,
        settings.template,
        settings.pressures,
        new_eos;
        kwargs...,
    )
end

function launchjob(calc::T, configfile::AbstractString; kwargs...) where {T<:ScfOrOptim}
    settings = load_settings(configfile)
    inputs =
        @. settings.dirs * '/' * (T <: SelfConsistentField ? "scf" : "vc-relax") * ".in"
    outputs = map(Base.Fix2(replace, ".in" => ".out"), inputs)
    return launchjob(calc, outputs, inputs, settings.manager.np, settings.bin; kwargs...)
end

function (step::Step{<:ScfOrOptim,Action{:fit_eos}})(
    outputs,
    trial_eos::EquationOfState,
    fit_e::Bool = true,
)
    results = map(outputs) do output
        _readdata(step, read(output, String))  # volume => energy
    end
    results = filter(Base.Fix2(!==, nothing), results)
    if length(results) <= 5
        @info "pressures <= 5 may give unreliable results, run more if possible!"
    end
    if fit_e
        return lsqfit(trial_eos(Energy()), first.(results), last.(results))
    else
        return lsqfit(trial_eos(Pressure()), first.(results), last.(results))
    end
end

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

"""
    postprocess(calc::Union{SelfConsistentField,VariableCellOptimization}, outputs, trial_eos::EquationOfState, fit_e::Bool = true)

Return the fitted equation of state from `outputs` and a `trial_eos`. Use `fit_e` to determine fit ``E(V)`` or ``P(V)``.
"""
function postprocess(
    calc::Union{SelfConsistentField,VariableCellOptimization},
    outputs,
    trial_eos::EquationOfState,
    fit_e::Bool = true,
)
    STEP_TRACKER[calc isa SelfConsistentField ? 3 : 6] =
        Context(nothing, outputs, Succeeded(), now(), Step(calc, ANALYSE_OUTPUT))
    return Step(calc, FIT_EOS)(outputs, trial_eos, fit_e)
end
"""
    postprocess(calc::Union{SelfConsistentField,VariableCellOptimization}, configfile)

Do the same thing of `postprocess`, but from a configuration file.
"""
function postprocess(calc::SelfConsistentField, configfile)
    settings = load_settings(configfile)
    inputs = settings.dirs .* "/scf.in"
    outputs = map(Base.Fix2(replace, ".in" => ".out"), inputs)
    return postprocess(calc, outputs, settings.trial_eos)
end
function postprocess(::VariableCellOptimization, configfile)
    settings = load_settings(configfile)
    inputs = settings.dirs .* "/vc-relax.in"
    outputs = map(Base.Fix2(replace, ".in" => ".out"), inputs)
    new_eos = postprocess(SelfConsistentField(), configfile)
    return postprocess(VariableCellOptimization(), outputs, new_eos)
end

function prep_potential(template)
    required = getpotentials(template)
    path = getpotentialdir(template)
    return map(required) do potential
        download_potential(potential, path)
    end
end

STEP_TRACKER = [
    Context(nothing, nothing, Pending(), now(), Step(SelfConsistentField(), PREPARE_INPUT)),
    Context(nothing, nothing, Pending(), now(), Step(SelfConsistentField(), LAUNCH_JOB)),
    Context(
        nothing,
        nothing,
        Pending(),
        now(),
        Step(SelfConsistentField(), ANALYSE_OUTPUT),
    ),
    Context(
        nothing,
        nothing,
        Pending(),
        now(),
        Step(VariableCellOptimization(), PREPARE_INPUT),
    ),
    Context(
        nothing,
        nothing,
        Pending(),
        now(),
        Step(VariableCellOptimization(), LAUNCH_JOB),
    ),
    Context(
        nothing,
        nothing,
        Pending(),
        now(),
        Step(VariableCellOptimization(), ANALYSE_OUTPUT),
    ),
]

function alert_pressures(pressures)
    if length(pressures) <= 5
        @info "pressures <= 5 may give unreliable results, consider more if possible!"
    end
    if minimum(pressures) >= zero(eltype(pressures))
        @warn "for better fitting, we need at least 1 negative pressure!"
    end
end # function alert_pressures

function _set_pressure_volume end

function _set_structure end

function getpotentials end

function getpotentialdir end

function download_potential end

function _readdata end

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

function load_settings(configfile)
    settings = load(configfile)
    _check_settings(settings)  # Errors will be thrown if exist
    return _expand_settings(settings)
end # function load_settings

include("QuantumESPRESSO.jl")

function Base.show(io::IO, eos::EquationOfState)  # Ref: https://github.com/mauro3/Parameters.jl/blob/3c1d72b/src/Parameters.jl#L542-L549
    if get(io, :compact, false)
        Base.show_default(IOContext(io, :limit => true), eos)
    else
        # just dumping seems to give ok output, in particular for big data-sets:
        T = typeof(eos)
        println(io, T)
        for f in fieldnames(T)
            println(io, " ", f, " = ", getfield(eos, f))
        end
    end
end # function Base.show

end
