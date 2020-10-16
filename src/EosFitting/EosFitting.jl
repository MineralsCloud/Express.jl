module EosFitting

using AbInitioSoftwareBase: loadfile
using AbInitioSoftwareBase.Inputs: Input, inputstring, writeinput
using Compat: isnothing
using EquationsOfStateOfSolids.Collections: EnergyEOS, PressureEOS
using EquationsOfStateOfSolids.Fitting: eosfit
using EquationsOfStateOfSolids.Volume: mustfindvolume
using OptionalArgChecks: @argcheck

using ..Express: ElectronicStructure, Optimization

import AbInitioSoftwareBase.Inputs: set_press_vol
import ..Express.Jobs: launchjob

export SelfConsistentField,
    StructuralOptimization,
    VariableCellOptimization,
    load_settings,
    set_press_vol,
    inputstring,
    prepare,
    finish,
    fiteos,
    writeinput,
    launchjob

struct SelfConsistentField <: ElectronicStructure end
struct StructuralOptimization <: Optimization end
struct VariableCellOptimization <: Optimization end

const ScfOrOptim = Union{SelfConsistentField,Optimization}

"""
    set_press_vol(template::Input, pressure, eos::EquationOfState; volume_scale = (0.5, 1.5))

Set the volume of `template` at a `pressure` according to `eos`.

The `volume_scale` gives a trial of the minimum and maximum scales for the `eos`. It
times the zero-pressure volume of the `eos` will be the trial volumes.
"""
function set_press_vol(
    template::Input,
    pressure,
    eos::PressureEOS;
    volume_scale = (0.5, 1.5),
)::Input
    volume = mustfindvolume(eos, pressure; volume_scale = volume_scale)
    return set_press_vol(template, pressure, volume)
end

"""
    prepare(calc, files, template::Input, pressures, trial_eos::EquationOfState; dry_run = false, kwargs...)
    prepare(calc, files, templates, pressures, trial_eos::EquationOfState; dry_run = false, kwargs...)

Prepare the input `files` from a certain `template` / a series of `templates` at `pressures` from a `trial_eos`.

Set `dry_run = true` to preview changes.
"""
function prepare(
    calc::ScfOrOptim,
    files,
    templates,
    pressures,
    trial_eos::PressureEOS;
    dry_run = false,
    kwargs...,
)
    alert_pressures(pressures)
    objects = map(files, templates, pressures) do file, template, pressure
        object = preset_template(calc, template)
        object = set_press_vol(object, pressure, trial_eos; kwargs...)
        writeinput(file, object, dry_run)
        object
    end
    return objects
end
prepare(
    calc::ScfOrOptim,
    files,
    template::Input,
    pressures,
    trial_eos::PressureEOS;
    kwargs...,
) = prepare(calc, files, fill(template, size(files)), pressures, trial_eos; kwargs...)
"""
    prepare(calc, configfile; kwargs...)

Do the same thing of `prepare`, but from a configuration file.
"""
function prepare(calc::SelfConsistentField, configfile; kwargs...)
    settings = load_settings(configfile)
    inputs = settings.dirs .* "/scf.in"
    return prepare(
        calc,
        inputs,
        settings.template,
        settings.pressures,
        PressureEOS(settings.trial_eos);
        kwargs...,
    )
end
function prepare(calc::VariableCellOptimization, configfile; kwargs...)
    settings = load_settings(configfile)
    inputs = settings.dirs .* "/vc-relax.in"
    new_eos = finish(SelfConsistentField(), configfile)
    return prepare(
        calc,
        inputs,
        settings.template,
        settings.pressures,
        PressureEOS(new_eos);
        kwargs...,
    )
end

function launchjob(::T, configfile; kwargs...) where {T<:ScfOrOptim}
    settings = load_settings(configfile)
    inputs =
        @. settings.dirs * '/' * (T <: SelfConsistentField ? "scf" : "vc-relax") * ".in"
    outputs = map(Base.Fix2(replace, ".in" => ".out"), inputs)
    return launchjob(outputs, inputs, settings.manager.np, settings.bin; kwargs...)
end

"""
    fiteos(calc, outputs, trial_eos::EquationOfState, fit_energy::Bool = true)

Fit an equation of state from `outputs` and a `trial_eos`. Use `fit_e` to determine fit ``E(V)`` or ``P(V)``.
"""
function fiteos(calc::ScfOrOptim, outputs, trial_eos::EnergyEOS)
    data =
        filter(!isnothing, [_readoutput(calc, read(output, String)) for output in outputs])  # volume => energy
    if length(data) <= 5
        @info "pressures <= 5 may give unreliable results, run more if possible!"
    end
    return eosfit(trial_eos, first.(data), last.(data))
end

"""
    finish(calc, outputs, trial_eos::EquationOfState, fit_energy::Bool = true)

Return the fitted equation of state from `outputs` and a `trial_eos`. Use `fit_e` to determine fit ``E(V)`` or ``P(V)``.
"""
function finish(calc::ScfOrOptim, outputs, trial_eos::EnergyEOS)
    return fiteos(calc, outputs, trial_eos)
end
"""
    finish(calc, configfile)

Do the same thing of `finish`, but from a configuration file.
"""
function finish(calc::SelfConsistentField, configfile)
    settings = load_settings(configfile)
    inputs = settings.dirs .* "/scf.in"
    outputs = map(Base.Fix2(replace, ".in" => ".out"), inputs)
    return finish(calc, outputs, EnergyEOS(settings.trial_eos))
end
function finish(::VariableCellOptimization, configfile)
    settings = load_settings(configfile)
    inputs = settings.dirs .* "/vc-relax.in"
    outputs = map(Base.Fix2(replace, ".in" => ".out"), inputs)
    new_eos = finish(SelfConsistentField(), configfile)
    return finish(VariableCellOptimization(), outputs, EnergyEOS(new_eos))
end

function alert_pressures(pressures)
    if length(pressures) <= 5
        @info "pressures <= 5 may give unreliable results, consider more if possible!"
    end
    if minimum(pressures) >= zero(eltype(pressures))
        @warn "for better fitting, we need at least 1 negative pressure!"
    end
end

function preset_template end

function _readoutput end

function _expand_settings end

function _check_software_settings end

function _check_settings(settings)
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
end

function load_settings(configfile)
    settings = loadfile(configfile)
    _check_settings(settings)  # Errors will be thrown if exist
    return _expand_settings(settings)
end

include("QuantumESPRESSO.jl")

end
