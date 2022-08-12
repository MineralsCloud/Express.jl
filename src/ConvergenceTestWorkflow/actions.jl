using AbInitioSoftwareBase: save, load, extension
using AbInitioSoftwareBase.Inputs: Input, writetxt, getpseudodir, getpotentials
using Dates: now, format
using ExpressBase: Action, calculation
using SimpleWorkflows: Job
using Unitful: ustrip, unit

using ..Shell: distprocs
using .Config: ExpandConfig

struct MakeInput{T} <: Action{T} end
function (x::MakeInput)(file, template::Input, args...)
    input = x(template, args...)
    mkpath(dirname(file))  # In case its parent directory is not created
    writetxt(file, input)
    return input
end

function buildjob(x::MakeInput{Scf}, cfgfile)
    dict = load(cfgfile)
    config = ExpandConfig{Scf}()(dict)
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

struct GetData{T} <: Action{T} end
function (x::GetData)(outputs)
    raw = (parseoutput(output) for output in outputs)  # `ntuple` cannot work with generators
    return collect(Iterators.filter(x -> x !== nothing, raw))  # A vector of pairs
end

function buildjob(x::GetData{T}, cfgfile) where {T}
    dict = load(cfgfile)
    config = ExpandConfig{T}()(dict)
    return Job(function ()
        data = x(last.(config.files))
        saved = Dict("results" => (ustrip ∘ last).(data))
        save(config.save_raw, saved)
        return data
    end)
end

struct TestConvergence{T} <: Action{T} end
(x::TestConvergence)(data) = isconvergent(data)

function buildjob(x::TestConvergence{T}, cfgfile) where {T}
    dict = load(cfgfile)
    config = ExpandConfig{T}()(dict)
    return Job(function ()
        data = GetData{T}()(last.(config.files))
        saved = Dict("results" => (ustrip ∘ last).(data))
        save(config.save_raw, saved)
        return x(data)
    end)
end
