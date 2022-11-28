module Config

using AbInitioSoftwareBase.Commands: CommandConfig
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
using ExpressBase: Action, Scf
using Unitful: Quantity

using ...Config: Directory, getfiles, _uparse, @sp

import Configurations: from_dict

@sp Pressures "GPa" "pressures" begin
    function (numbers, unit)
        if length(numbers) <= 5
            @info "less than 6 pressures may not fit accurately, consider adding more!"
        end
        if minimum(numbers * unit) >= 0 * unit  # `numbers` may have eltype `Any`
            @warn "for better fitting result, provide at least 1 negative pressure!"
        end
    end
end

@sp Volumes "bohr^3" "volumes" begin
    function (numbers, _)
        if length(numbers) <= 5
            @info "less than 6 volumes may not fit accurately, consider adding more!"
        end
    end
end

@option struct TrialEquationOfState
    type::String
    params::Vector{Quantity{Float64}}
end

@option struct Save
    ev::String = "ev.json"
    eos::String = "eos.jld2"
    wf::String = "wf.jld2"
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

struct ExpandConfig{T} <: Action{T} end
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
function (::ExpandConfig{T})(dir::Directory, fixed::Union{Pressures,Volumes}) where {T}
    return map(fixed.numbers) do number
        getfiles(dir, number, string(nameof(T)))
    end
end
function (::ExpandConfig)(save::Save)
    keys = fieldnames(Save)
    values = (abspath(expanduser(getfield(save, key))) for key in keys)
    return (; zip(keys, values)...)
end
function (x::ExpandConfig)(config::RuntimeConfig)
    return (
        template=x(config.template),
        trial_eos=PressureEquation(x(config.trial_eos)),
        fixed=x(config.fixed),
        files=x(config.dir, config.fixed),
        scffiles=ExpandConfig{Scf}()(config.dir, config.fixed),
        save=x(config.save),
        cli=config.cli,
    )
end

function from_dict(
    ::Type{TrialEquationOfState}, ::OptionField{:params}, ::Type{<:Quantity}, param
)
    return eval(_uparse(string(param)))
end

end
