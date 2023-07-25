using AbInitioSoftwareBase: Input, writetxt
using ExpressBase: Action, SCF
using ExpressBase.Files: save, load
using EasyJobsBase: Job
using Unitful: ustrip

using ..Express: distribute_procs
using .Config: ExpandConfig

import ..Express: jobify

struct CreateInput{T} <: Action{T}
    calculation::T
end

function jobify(x::CreateInput{SCF}, cfgfile)
    dict = load(cfgfile)
    config = ExpandConfig{SCF}()(dict)
    inputs = first.(config.files)
    if config.parameters isa AbstractVector{<:Tuple}
        return map(inputs, config.parameters) do input, (mesh, shift)
            Job(() -> x(input, config.template, mesh, shift, "Y-m-d_H:M:S"))
        end
    else
        return map(inputs, config.parameters) do input, energy
            Job(() -> x(input, config.template, energy, "Y-m-d_H:M:S"))
        end
    end
end

function parseoutput end

struct ExtractData{T} <: Action{T}
    calculation::T
end

function jobify(x::ExtractData{T}, cfgfile) where {T}
    dict = load(cfgfile)
    config = ExpandConfig{T}()(dict)
    return Job(function ()
        data = x(last.(config.files))
        saved = Dict("results" => (ustrip ∘ last).(data))
        save(config.save_raw, saved)
        return data
    end)
end

struct TestConvergence{T} <: Action{T}
    calculation::T
end
(x::TestConvergence)(data) = isconvergent(data)

function jobify(x::TestConvergence{T}, cfgfile) where {T}
    dict = load(cfgfile)
    config = ExpandConfig{T}()(dict)
    return Job(function ()
        data = ExtractData{T}()(last.(config.files))
        saved = Dict("results" => (ustrip ∘ last).(data))
        save(config.save_raw, saved)
        return x(data)
    end)
end

function isconvergent(a::AbstractVector)
    terms = abs.(diff(a))
    x, y, z = last(terms, 3)
    return all(0 <= r < 1 for r in (y / x, z / y))
end
