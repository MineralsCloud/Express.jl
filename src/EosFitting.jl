module EosFitting

using AbInitioSoftwareBase: loadfile
using AbInitioSoftwareBase.Inputs: Input, inputstring, writeinput
using Compat: isnothing
using EquationsOfStateOfSolids.Collections: EquationOfStateOfSolids, EnergyEOS, PressureEOS
using EquationsOfStateOfSolids.Fitting: eosfit
using EquationsOfStateOfSolids.Volume: mustfindvolume
using OptionalArgChecks: @argcheck

using ..Express: ElectronicStructure, Optimization

import AbInitioSoftwareBase.Inputs: set_press_vol

export SelfConsistentField,
    StructuralOptimization,
    VariableCellOptimization,
    load_settings,
    set_press_vol,
    inputstring,
    prepare_input,
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
function set_press_vol(template::Input, pressure, eos::PressureEOS)::Input
    volume = mustfindvolume(eos, pressure; volume_scale = vscaling())
    return set_press_vol(template, pressure, volume)
end

"""
    prepare(calc, files, template::Input, pressures, trial_eos::EquationOfState; dry_run = false, kwargs...)
    prepare(calc, files, templates, pressures, trial_eos::EquationOfState; dry_run = false, kwargs...)

Prepare the input `files` from a certain `template` / a series of `templates` at `pressures` from a `trial_eos`.

Set `dry_run = true` to preview changes.
"""
function prepare_input(calc::ScfOrOptim)
    function _prepare_input(file, template::Input, pressure, eos_or_volume; kwargs...)
        object = customize(standardize(template, calc), pressure, eos_or_volume; kwargs...)
        writeinput(file, object)
        return object
    end
    function _prepare_input(files, templates, pressures, eos_or_volumes; kwargs...)
        _alert(pressures)
        if templates isa Input
            templates = fill(templates, size(files))
        end
        objects = if eos_or_volumes isa EquationOfStateOfSolids
            map(files, templates, pressures) do file, template, pressure
                _prepare_input(file, template, pressure, eos_or_volumes; kwargs...)
            end
        else
            map(files, templates, pressures, eos_or_volumes) do file, template, pressure, volume
                _prepare_input(file, template, pressure, volume; kwargs...)
            end
        end
        return objects
    end
    function _prepare_input(cfgfile; kwargs...)
        settings = load_settings(cfgfile)
        inputs = settings.dirs .* settings.name
        _param(::SelfConsistentField) = settings.trial_eos
        _param(::VariableCellOptimization) = finish(SelfConsistentField(), cfgfile)
        return _prepare_input(
            inputs,
            settings.template,
            settings.pressures,
            PressureEOS(_param(calc));
            kwargs...,
        )
    end
    return _prepare_input
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

function _alert(pressures)
    if length(pressures) <= 5
        @info "pressures <= 5 may give unreliable results, consider more if possible!"
    end
    if minimum(pressures) >= zero(eltype(pressures))
        @warn "for better fitting, we need at least 1 negative pressure!"
    end
end

function standardize end

function customize end

function _readoutput end

function _expand_settings end

function _check_software_settings end

vscaling()::NTuple{2,<:AbstractFloat} = (0.5, 1.5)

function _check_settings(settings)
    map(("template", "pressures", "trial_eos", "dir")) do key
        @argcheck haskey(settings, key)
    end
    _check_software_settings(settings["qe"])
    @argcheck isdir(settings["dir"])
    @argcheck isfile(settings["template"])
    _alert(settings["pressures"])
    map(("type", "parameters", "units")) do key
        @argcheck haskey(settings["trial_eos"], key)
    end
end

function load_settings(configfile)
    settings = loadfile(configfile)
    _check_settings(settings)  # Errors will be thrown if exist
    return _expand_settings(settings)
end

end
