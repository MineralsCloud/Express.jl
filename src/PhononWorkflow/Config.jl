module Config

using Configurations: from_dict, @option
using EasyConfig: Config as Conf
using ExpressBase: Calculation
using ExpressBase.Config: AbstractConfig, SoftwareConfig, SamplingPoints, IO, list_io
using Unitful: FreeUnits

@option "pressures" struct Pressures <: SamplingPoints
    numbers::Vector{Float64}
    unit::FreeUnits
    Pressures(numbers, unit="GPa") = new(numbers, unit)
end

@option "volumes" struct Volumes <: SamplingPoints
    numbers::Vector{Float64}
    unit::FreeUnits
    Volumes(numbers, unit="bohr^3") = new(numbers, unit)
end

@option struct Data <: AbstractConfig
    static::String = "energies.json"
end

@option struct StaticConfig <: AbstractConfig
    recipe::String
    templates::Vector{String}
    at::Union{Pressures,Volumes}
    io::IO = IO()
    data::Data = Data()
    cli::SoftwareConfig
    function StaticConfig(recipe, templates, at, io, data, cli)
        @assert recipe in ("phonon dispersion", "vdos")
        for template in templates
            if !isfile(template)
                @warn "I cannot find template file `$template`!"
            end
        end
        return new(recipe, templates, at, io, data, cli)
    end
end

function _update!(conf::Conf, at::Union{Pressures,Volumes})
    conf.at = collect(number for number in at)
    return conf
end
function _update!(conf::Conf, io::IO, at::Union{Pressures,Volumes})
    conf.io = collect(
        list_io(io, number, string(nameof(typeof(conf.calculation)))) for
        number in at.numbers
    )
    return conf
end
function _update!(conf::Conf, data::Data)
    conf.data.static = abspath(expanduser(data.static))
    return conf
end

function expand(config::StaticConfig, calculation::Calculation)
    conf = Conf()
    conf.cli = config.cli
    conf.calculation = calculation
    _update!(conf, config.templates)
    _update!(conf, config.at)
    _update!(conf, config.io, config.at)
    _update!(conf, config.data)
    return conf
end

end
