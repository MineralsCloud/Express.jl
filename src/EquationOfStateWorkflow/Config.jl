module Config

using Configurations: OptionField, @option
using EasyConfig: Config as Conf
using EquationsOfStateOfSolids:
    Murnaghan1st,
    Murnaghan2nd,
    BirchMurnaghan2nd,
    BirchMurnaghan3rd,
    BirchMurnaghan4th,
    PoirierTarantola2nd,
    PoirierTarantola3rd,
    PoirierTarantola4th,
    Vinet,
    PressureEquation
using ExpressBase: Action, Calculation, SelfConsistentField
using ExpressBase.Config:
    AbstractConfig,
    SoftwareConfig,
    SamplingPoints,
    IO,
    Subdirectory,
    InputFile,
    OutputFile,
    list_io,
    _uparse
using Unitful: Quantity, FreeUnits
using UnitfulParsableString  # Override `string`

import Configurations: from_dict

@option "pressures" struct Pressures <: SamplingPoints
    numbers::Vector{Float64}
    unit::FreeUnits
    function Pressures(numbers, unit="GPa")
        if length(numbers) <= 5
            @info "less than 6 pressures may not fit accurately, consider adding more!"
        end
        if minimum(numbers * unit) >= 0 * unit  # `numbers` may have eltype `Any`
            @warn "for better fitting result, provide at least 1 negative pressure!"
        end
        return new(numbers, unit)
    end
end

@option "volumes" struct Volumes <: SamplingPoints
    numbers::Vector{Float64}
    unit::FreeUnits
    function Volumes(numbers, unit="bohr^3")
        if length(numbers) <= 5
            @info "less than 6 volumes may not fit accurately, consider adding more!"
        end
        return new(numbers, unit)
    end
end

@option struct TrialEquationOfState <: AbstractConfig
    type::String
    params::Vector{Quantity{Float64}}
end

@option struct Data <: AbstractConfig
    raw::String = "volumes_energies.json"
    eos_params::String = "eos_params.json"
end

@option struct StaticConfig <: AbstractConfig
    recipe::String
    template::String
    trial_eos::TrialEquationOfState
    at::Union{Pressures,Volumes}
    io::IO = IO()
    data::Data = Data()
    cli::CommandConfig
    function StaticConfig(recipe, template, trial_eos, at, io, save, cli)
        @assert recipe in ("eos",)
        if !isfile(template)
            @warn "I cannot find template file `$template`!"
        end
        return new(recipe, template, trial_eos, at, io, save, cli)
    end
end

function expand!(conf::Conf, trial_eos::TrialEquationOfState, ::Calculation)
    type = filter(c -> isletter(c) || isdigit(c), lowercase(trial_eos.type))
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
    conf.trial_eos = T(trial_eos.params...)
    return conf
end
function expand!(conf::Conf, io::IO, at::Union{Pressures,Volumes}, ::Calculation)
    conf.io = collect(
        list_io(io, number, string(nameof(typeof(conf.calculation)))) for
        number in at.numbers
    )
    return conf
end
function expand!(conf::Conf, data::Data, ::Calculation)
    conf.data.raw = abspath(expanduser(data.raw))
    conf.data.eos_params = abspath(expanduser(data.eos_params))
    return conf
end

function expand(config::StaticConfig, calculation::Calculation)
    conf = Conf()
    conf.cli = config.cli
    expand!(conf, config.template, calculation)
    expand!(conf, config.trial_eos, calculation)
    expand!(conf, config.io, config.at, calculation)
    expand!(conf, config.data, calculation)
    return conf
end

from_dict(::Type{TrialEquationOfState}, ::OptionField{:params}, ::Type{<:Quantity}, param) =
    eval(_uparse(string(param)))

end
