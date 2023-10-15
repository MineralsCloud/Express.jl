module Config

using Configurations: OptionField, @option
using EasyConfig: Config as Conf
using ExpressBase: Calculation
using ExpressBase.Config:
    AbstractConfig, SoftwareConfig, SamplingPoints, IO, list_io, _uparse
using Unitful: FreeUnits, Quantity

import Configurations: from_dict

@option "ecut" struct CutoffEnergies <: SamplingPoints
    numbers::Vector{Float64}
    unit::FreeUnits
    CutoffEnergies(numbers, unit="Ry") = new(numbers, unit)
end

@option "kmesh" struct MonkhorstPackGrids <: AbstractConfig
    meshes::AbstractVector{<:AbstractVector{<:Integer}}
    shifts::AbstractVector{<:AbstractVector{Bool}} = fill(falses(3), length(meshes))
    function MonkhorstPackGrids(meshes, shifts)
        if length(meshes) != length(shifts)
            throw(DimensionMismatch("`meshes` and `shifts` should have the same length!"))
        end
        for (mesh, shift) in zip(meshes, shifts)
            @assert all(mesh .>= 1)
            @assert all(0 .<= shift .<= 1)
        end
        return new(meshes, shifts)
    end
end

@option struct Data <: AbstractConfig
    raw::String = "energies.json"
end

@option struct StaticConfig <: AbstractConfig
    recipe::String
    template::String
    with::Union{CutoffEnergies,MonkhorstPackGrids}
    threshold::Quantity
    io::IO = IO()
    data::Data = Data()
    cli::SoftwareConfig
    function StaticConfig(recipe, template, with, threshold, io, data, cli)
        @assert recipe in ("ecut", "kmesh")
        if !isfile(template)
            @warn "I cannot find template file `$template`!"
        end
        return new(recipe, template, with, threshold, io, data, cli)
    end
end

function _update!(conf::Conf, with::CutoffEnergies)
    conf.with = collect(number for number in with)
    return conf
end
function _update!(conf::Conf, with::MonkhorstPackGrids)
    conf.with = collect((mesh, shift) for (mesh, shift) in zip(with.meshes, with.shifts))
    return conf
end
function _update!(conf::Conf, io::IO, energies::CutoffEnergies)
    conf.io = collect(
        list_io(io, number, string(nameof(typeof(conf.calculation)))) for
        number in energies.numbers
    )
    return conf
end
function _update!(conf::Conf, io::IO, grids::MonkhorstPackGrids)
    conf.io = collect(
        list_io(
            io,
            join(mesh, '×') * '_' * join(shift, '+'),
            string(nameof(typeof(conf.calculation))),
        ) for (mesh, shift) in zip(grids.meshes, grids.shifts)
    )
    return conf
end
function _update!(conf::Conf, data::Data)
    conf.data.raw = abspath(expanduser(data.raw))
    return conf
end
function _update!(conf::Conf, threshold::Quantity)
    conf.threshold = threshold
    return conf
end

function expand(config::StaticConfig, calculation::Calculation)
    conf = Conf()
    conf.cli = config.cli
    conf.calculation = calculation
    _update!(conf, config.template)
    _update!(conf, config.with)
    _update!(conf, config.io, config.with)
    _update!(conf, config.data)
    _update!(conf, config.threshold)
    return conf
end

from_dict(
    ::Type{StaticConfig}, ::OptionField{:criteria}, ::Type{Quantity}, str::AbstractString
) = eval(_uparse(str))

end
