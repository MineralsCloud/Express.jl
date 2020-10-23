module EosFitting

using AbInitioSoftwareBase: loadfile
using AbInitioSoftwareBase.CLI: MpiExec
using AbInitioSoftwareBase.Inputs: Input, inputstring, writeinput
using Compat: isnothing
using EquationsOfStateOfSolids.Collections:
    EquationOfStateOfSolids,
    EnergyEOS,
    PressureEOS,
    Murnaghan,
    BirchMurnaghan2nd,
    BirchMurnaghan3rd,
    BirchMurnaghan4th,
    Vinet
using EquationsOfStateOfSolids.Volume: mustfindvolume
using Mustache: render
using SimpleWorkflow: ExternalAtomicJob, InternalAtomicJob, Script, chain
using Unitful: uparse
import Unitful
import UnitfulAtomic
using UrlDownload: File, URL, urldownload

using ..Express: ElectronicStructure, Optimization

import EquationsOfStateOfSolids.Fitting: eosfit
import AbInitioSoftwareBase.Inputs: set_press_vol

export SelfConsistentField,
    StructureOptimization,
    VariableCellOptimization,
    load_settings,
    inputstring,
    makeinput,
    readoutput,
    eosfit,
    writeinput

struct SelfConsistentField <: ElectronicStructure end
struct StructureOptimization <: Optimization end
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
function makeinput(calc::ScfOrOptim)
    function _makeinput(file, template::Input, pressure, eos_or_volume; kwargs...)
        object = customize(standardize(template, calc), pressure, eos_or_volume; kwargs...)
        writeinput(file, object)
        return object
    end
    function _makeinput(files, templates, pressures, eos_or_volumes; kwargs...)
        _alert(pressures)
        if templates isa Input
            templates = fill(templates, size(files))
        end
        objects = if eos_or_volumes isa EquationOfStateOfSolids
            map(files, templates, pressures) do file, template, pressure
                _makeinput(file, template, pressure, eos_or_volumes; kwargs...)
            end
        else
            map(files, templates, pressures, eos_or_volumes) do file, template, pressure, volume
                _makeinput(file, template, pressure, volume; kwargs...)
            end
        end
        return objects
    end
    function _makeinput(cfgfile; kwargs...)
        settings = load_settings(cfgfile)
        inputs = joinpath.(settings.dirs, settings.name)
        eos = PressureEOS(
            calc isa SelfConsistentField ? settings.trial_eos :
            eosfit(SelfConsistentField())(cfgfile),
        )
        return _makeinput(inputs, settings.templates, settings.pressures, eos; kwargs...)
    end
end

abstract type JobPackaging end
struct JobOfTasks <: JobPackaging end
struct ArrayOfJobs <: JobPackaging end

function buildjob(::typeof(makeinput), calc::ScfOrOptim)
    function _buildjob(file, template::Input, pressure, eos_or_volume; kwargs...)
        f = makeinput(calc)
        return InternalAtomicJob(
            () -> f(file, template, pressure, eos_or_volume; kwargs...),
            "Prepare $calc input for pressure $pressure",
        )
    end
    function _buildjob(files, templates, pressures, eos_or_volumes, ::JobOfTasks; kwargs...)
        f = makeinput(calc)
        return InternalAtomicJob(
            () -> f(files, templates, pressures, eos_or_volumes; kwargs...),
            "Prepare $calc inputs for pressures $pressures",
        )
    end
    function _buildjob(
        files,
        templates,
        pressures,
        eos_or_volumes,
        ::ArrayOfJobs;
        kwargs...,
    )
        if templates isa Input
            templates = fill(templates, size(files))
        end
        if eos_or_volumes isa EquationOfStateOfSolids
            map(files, templates, pressures) do file, template, pressure
                _buildjob(file, template, pressure, eos_or_volumes; kwargs...)
            end
        else
            map(
                files,
                templates,
                pressures,
                eos_or_volumes,
            ) do file, template, pressure, volume
                _buildjob(file, template, pressure, volume; kwargs...)
            end
        end
    end
end
function buildjob(::typeof(eosfit), calc::ScfOrOptim)
    function _buildjob(args...)
        return InternalAtomicJob(
            () -> eosfit(calc)(args...),
            "Prepare $calc inputs for pressures ",
        )
    end
end

function readoutput(calc::ScfOrOptim)
    function _readoutput(str::AbstractString, parser = nothing)
        if isnothing(parser)
            return str
        else
            return parser(str, calc)  # `parseoutput` will be used here
        end
    end
    function _readoutput(url_or_file::Union{URL,File}, parser = nothing)
        str = urldownload(url_or_file, true; parser = String)
        return _readoutput(str, parser)
    end
    function _readoutput(file, parser = nothing)
        open(file, "r") do io
            str = read(io, String)
            return _readoutput(str, parser)
        end
    end
end

function makescript(template, view)
    map((:press, :nprocs, :in, :out, :script)) do key
        @assert haskey(view, key)
    end
    str = render(template, view)
    return Script(str, view[:script])
end
makescript(template, args::Pair...) = makescript(template, Dict(args))
makescript(template; kwargs...) = makescript(template, Dict(kwargs))

function distprocs(nprocs, njobs)
    quotient, remainder = divrem(nprocs, njobs)
    if !iszero(remainder)
        @warn "The processes are not fully balanced! Consider the number of subjobs!"
    end
    return quotient
end

function buildjob(::ScfOrOptim)
    function _buildjob(outputs, inputs, np, exe; kwargs...)
        # `map` guarantees they are of the same size, no need to check.
        n = distprocs(np, length(inputs))
        subjobs = map(outputs, inputs) do output, input
            f = MpiExec(np; kwargs...) ∘ exe
            cmd = f(stdin = input, stdout = output)
            ExternalAtomicJob(cmd)
        end
        return subjobs
    end
    function _buildjob(template, view)
        ExternalAtomicJob(makescript(template, view))
    end
end

"""
    fiteos(calc, outputs, trial_eos::EquationOfState, fit_energy::Bool = true)

Fit an equation of state from `outputs` and a `trial_eos`. Use `fit_e` to determine fit ``E(V)`` or ``P(V)``.
"""
function eosfit(calc::ScfOrOptim)
    function _eosfit(outputs, trial_eos::EnergyEOS)
        reader = readoutput(calc)
        raw = (reader(output, parseoutput) for output in outputs)  # `ntuple` cannot work with generators
        data = collect(Iterators.filter(!isnothing, raw))  # A vector of pairs
        if length(data) <= 5
            @info "pressures <= 5 may give unreliable results, run more if possible!"
        end
        return eosfit(trial_eos, first.(data), last.(data))
    end
    function _eosfit(cfgfile)
        settings = load_settings(cfgfile)
        inputs = settings.dirs .* settings.name
        outputs = map(Base.Fix2(replace, ".in" => ".out"), inputs)
        eos = EnergyEOS(
            calc isa SelfConsistentField ? settings.trial_eos :
            eosfit(SelfConsistentField())(cfgfile),
        )
        return _eosfit(outputs, eos)
    end
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

function parseoutput end

function expand_settings end

function check_software_settings end

vscaling()::NTuple{2,<:AbstractFloat} = (0.5, 1.5)

const UNIT_CONTEXT = [Unitful, UnitfulAtomic]

function expandeos(settings)
    type = string(lowercase(settings["type"]))
    T = if type in ("m", "murnaghan")
        Murnaghan
    elseif type in ("bm2", "birchmurnaghan2nd", "birch-murnaghan-2")
        BirchMurnaghan2nd
    elseif type in ("bm3", "birchmurnaghan3rd", "birch-murnaghan-3")
        BirchMurnaghan3rd
    elseif type in ("bm4", "birchmurnaghan4th", "birch-murnaghan-4")
        BirchMurnaghan4th
    elseif type in ("v", "vinet")
        Vinet
    end
    p = map(settings["parameters"]) do x
        @assert length(x) == 2
        first(x) * uparse(last(x); unit_context = UNIT_CONTEXT)
    end
    return T(p...)
end

function _check_settings(settings)
    map(("templates", "pressures", "trial_eos", "workdir")) do key
        @assert haskey(settings, key)
    end
    @assert isdir(settings["workdir"])
    @assert all(map(isfile, settings["templates"]))
    _alert(settings["pressures"]["values"])
    map(("type", "parameters")) do key
        @assert haskey(settings["trial_eos"], key)
    end
end

function load_settings(configfile)
    settings = loadfile(configfile)
    _check_settings(settings)  # Errors will be thrown if exist
    return expand_settings(settings)
end

end
