module Config

using Configurations: OptionField, @option
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
using ExpressBase: Action, SCF, CommandConfig
using ExpressBase.Config: SamplingPoints, IO, Subdirectory, list_io, _uparse
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

@option struct TrialEquationOfState
    type::String
    params::Vector{Quantity{Float64}}
end

@option struct Save
    raw_data::String = "volumes_energies.json"
    parameters::String = "parameters.json"
end

@option struct RuntimeConfig
    recipe::String
    template::String
    trial_eos::TrialEquationOfState
    fixed::Union{Pressures,Volumes}
    dir::Directory = Directory()
    save::Save = Save()
    cli::CommandConfig
    function RuntimeConfig(recipe, template, trial_eos, fixed, dir, save, cli)
        @assert recipe in ("eos",)
        if !isfile(template)
            @warn "I cannot find template file `$template`!"
        end
        return new(recipe, template, trial_eos, fixed, dir, save, cli)
    end
end

struct ExpandConfig{T} <: Action{T}
    calculation::T
end
function (::ExpandConfig)(trial_eos::TrialEquationOfState)
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
    return T(trial_eos.params...)
end
(::ExpandConfig)(data::Union{Pressures,Volumes}) = collect(datum for datum in data)
function (obj::ExpandConfig)(dir::Directory, fixed::Union{Pressures,Volumes})
    return map(fixed.numbers) do number
        getfiles(dir, number, string(obj.calculation))
    end
end
function (::ExpandConfig)(save::Save)
    keys = fieldnames(Save)
    values = (abspath(expanduser(getfield(save, key))) for key in keys)
    return (; zip(keys, values)...)
end
function (obj::ExpandConfig)(config::RuntimeConfig)
    return (
        template=obj(config.template),
        trial_eos=PressureEquation(obj(config.trial_eos)),
        fixed=obj(config.fixed),
        files=obj(config.dir, config.fixed),
        scffiles=ExpandConfig(SCF())(config.dir, config.fixed),
        save=obj(config.save),
        cli=config.cli,
    )
end

function from_dict(
    ::Type{TrialEquationOfState}, ::OptionField{:params}, ::Type{<:Quantity}, param
)
    return eval(_uparse(string(param)))
end

end
