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
    getparam
using EquationsOfStateOfSolids.Fitting: eosfit
using ExpressBase: Calculation, Action, Scf, Optimization, calculation
using ExpressBase.Files: save, load, extension
using Unitful: Pressure, Volume, ustrip, unit
using UnitfulParsableString: string

using .Config: ExpandConfig, Pressures, Volumes, _uparse

struct MakeInput{T} <: Action{T} end
(makeinput::MakeInput)(
    path, template::Input, pressure::Pressure, eos::PressureEquation, args...
) = makeinput(path, makeinput(template, pressure, eos, args...))
(makeinput::MakeInput)(path, template::Input, volume::Volume, args...) =
    makeinput(path, makeinput(template, volume, args...))
function (makeinput::MakeInput)(path, input::Input)
    if isfile(path)
        @warn "File $path already exists! It will be overwritten!"
    end
    mkpath(dirname(path))  # In case its parent directory is not created
    open(path, "w") do io
        print(io, input)
    end
    return input
end

struct ExtractData{T} <: Action{T} end

struct SaveVolumeEnergy{T} <: Action{T} end
function (::SaveVolumeEnergy{T})(path, data) where {T}
    data = sort(collect(data))  # In case the data is not sorted
    dict = Dict("volume" => first.(data), "energy" => last.(data))
    return save(path, dict)
end

struct FitEquationOfState{T} <: Action{T} end
function (fit::FitEquationOfState{T})(trial_eos::EnergyEquation) where {T}
    return function (data)
        data = sort(collect(data))  # In case the data is not sorted
        if length(data) <= 5
            @info "pressures <= 5 may give unreliable results, run more if possible!"
        end
        return eosfit(trial_eos, collect(first.(data)), collect(last.(data)))
    end
end

struct SaveParameters{T} <: Action{T} end
function (::SaveParameters)(path, parameters::Parameters)
    dict = Dict(
        "type" => string(typeof(parameters)),
        "params" => collect(
            string(getproperty(parameters, name)) for name in propertynames(parameters)
        ),
    )
    return save(path, dict)
end
(x::SaveParameters)(path, eos::EquationOfStateOfSolids) = x(path, getparam(eos))

struct LoadParameters{T} <: Action{T} end
function (::LoadParameters)(path)
    data = load(path)
    type, params = data["type"], data["params"]
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
        error("unsupported eos name `\"$type\"`!")
    end
    return T(map(_uparse, params)...)
end

include("think.jl")
