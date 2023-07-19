using AbInitioSoftwareBase: Input
using EquationsOfStateOfSolids:
    EquationOfStateOfSolids,
    EnergyEquation,
    PressureEquation,
    Parameters,
    Murnaghan1st,
    Murnaghan2nd,
    BirchMurnaghan2nd,
    BirchMurnaghan3rd,
    BirchMurnaghan4th,
    PoirierTarantola2nd,
    PoirierTarantola3rd,
    PoirierTarantola4th,
    Vinet,
    getparam,
    vsolve
using EquationsOfStateOfSolids.Fitting: eosfit
using ExpressBase:
    Calculation, Action, SCF, Optimization, DownloadPotentials, RunCmd, WriteInput
using ExpressBase.Config: _uparse
using ExpressBase.Files: save, load, extension
using Unitful: Pressure, Volume
using UnitfulParsableString: string

using .Config: Pressures, Volumes

struct ComputeVolume{T} <: Action{T}
    calculation::T
end
(obj::ComputeVolume)(pressure::Pressure, parameters::Parameters) =
    obj(pressure, PressureEquation(parameters))
function (obj::ComputeVolume)(pressure::Pressure, eos::PressureEquation)
    possible_volumes = vsolve(eos, pressure)
    return if length(possible_volumes) > 1
        _choose(possible_volumes, pressure, eos)
    else
        only(possible_volumes)
    end
end

struct CreateInput{T} <: Action{T}
    calculation::T
    template::Input
end

struct ExtractData{T} <: Action{T}
    calculation::T
end

struct SaveData{T} <: Action{T}
    calculation::T
    path::String
end
function (obj::SaveData)(data)
    data = sort(collect(data))  # In case the data is not sorted
    dict = Dict("volume" => (string ∘ first).(data), "energy" => (string ∘ last).(data))
    return save(obj.path, dict)
end

struct FitEquationOfState{T} <: Action{T}
    calculation::T
end
function (fit::FitEquationOfState)(trial_eos::EnergyEquation)
    return function (data)
        data = sort(collect(data))  # In case the data is not sorted
        if length(data) < length(fieldnames(typeof(trial_eos.param)))
            error("not enough data points to fit an EOS!")
        end
        if length(data) <= 5
            @info "pressures <= 5 may give unreliable results, run more if possible!"
        end
        return eosfit(trial_eos, first.(data), last.(data))
    end
end

struct SaveParameters{T} <: Action{T}
    calculation::T
    path::String
end
function (obj::SaveParameters)(parameters::Parameters)
    dict = Dict(
        "type" => if parameters isa Murnaghan1st
            "murnaghan"
        elseif parameters isa Murnaghan2nd
            "murnaghan2"
        elseif parameters isa BirchMurnaghan2nd
            "birchmurnaghan2"
        elseif parameters isa BirchMurnaghan3rd
            "birchmurnaghan3"
        elseif parameters isa BirchMurnaghan4th
            "birchmurnaghan4"
        elseif parameters isa PoirierTarantola2nd
            "poiriertarantola2"
        elseif parameters isa PoirierTarantola3rd
            "poiriertarantola3"
        elseif parameters isa PoirierTarantola4th
            "poiriertarantola4"
        elseif parameters isa Vinet
            "vinet"
        else
            error("unsupported EOS type `\"$parameters\"`!")
        end,
        "params" => collect(
            string(getproperty(parameters, name)) for name in propertynames(parameters)
        ),
    )
    return save(obj.path, dict)
end
(obj::SaveParameters)(eos::EquationOfStateOfSolids) = obj(getparam(eos))

function loadparameters(path)
    data = load(path)
    type, params = lowercase(data["type"]), data["params"]
    T = if type in ("m", "murnaghan")
        Murnaghan1st
    elseif type == "m2" || occursin("murnaghan2", type)
        Murnaghan2nd
    elseif type == "bm2" || occursin("birchmurnaghan2", type)
        BirchMurnaghan2nd
    elseif type == "bm3" || occursin("birchmurnaghan3", type)
        BirchMurnaghan3rd
    elseif type == "bm4" || occursin("birchmurnaghan4", type)
        BirchMurnaghan4th
    elseif type == "pt2" || occursin("poiriertarantola2", type)
        PoirierTarantola2nd
    elseif type == "pt3" || occursin("poiriertarantola3", type)
        PoirierTarantola3rd
    elseif type == "pt4" || occursin("poiriertarantola4", type)
        PoirierTarantola4th
    elseif type == "v" || occursin("vinet", type)
        Vinet
    else
        error("unsupported EOS type `\"$type\"`!")
    end
    return T(map(eval ∘ _uparse, params)...)
end

function _choose(possible_volumes, pressure, eos)
    v0 = getparam(eos).v0
    filtered = if pressure >= zero(pressure)  # If pressure is greater than zero,
        filter(<=(v0), possible_volumes)  # the volume could only be smaller than `v0`.
    else
        filter(v -> 1 < v / v0 <= 3, possible_volumes)
    end
    return only(filtered)
end
