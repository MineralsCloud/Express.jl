module Config

using Configurations: @option
using ExpressBase: Action
using Unitful: FreeUnits

using ...Config: SamplingPoints, Directory, getfiles, _uparse

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

@option struct RuntimeConfig
    recipe::String
    template::String
    fixed::Union{Pressures,Volumes}
    dir::Directory = Directory()
    save::Save = Save()
    cli::SoftwareConfig
    function RuntimeConfig(recipe, template, fixed, dir, save, cli)
        @assert recipe in ("static elasticity",)
        if !isfile(template)
            @warn "I cannot find template file `$template`!"
        end
        return new(recipe, template, fixed, dir, save, cli)
    end
end

struct ExpandConfig{T} <: Action{T} end
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

end
