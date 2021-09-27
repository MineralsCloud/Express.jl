module Recipes

using AbInitioSoftwareBase: load
using Serialization: deserialize
using SimpleWorkflows: AtomicJob, Workflow, run!, ▷, ⋲, ⋺

using ..Config: ExpandConfig
using ..EquationOfStateWorkflow:
    Scf, VcOptim, LogMsg, MakeInput, RunCmd, GetData, FitEos, buildjob

export build, run!

struct Recipe
    type::String
    config::Dict{String,Any}
end

function build(recipe::Recipe)
    type = Symbol(lowercase(recipe.config))
    return _build(recipe.config, Val(type))
end
function _build(config::AbstractDict, ::Val{:eos})
    if isfile(config["save"]["status"])
        w = deserialize(config["save"]["status"])
        typeassert(w, Workflow)
        return w
    else
        a = AtomicJob(() -> LogMsg{Scf}()(; start = true))
        b = buildjob(MakeInput{Scf}(), config)
        c = buildjob(RunCmd{Scf}(), config)
        d = buildjob(FitEos{Scf}(), config)
        f = AtomicJob(() -> LogMsg{Scf}()(; start = false))
        g = AtomicJob(() -> LogMsg{VcOptim}()(; start = true))
        h = buildjob(MakeInput{VcOptim}(), config)
        i = buildjob(RunCmd{VcOptim}(), config)
        j = buildjob(FitEos{VcOptim}(), config)
        l = AtomicJob(() -> LogMsg{VcOptim}()(; start = false))
        (((((((a ⋲ b) ▷ c) ⋺ d) ▷ f) ▷ g ⋲ h) ▷ i) ⋺ j) ▷ l
        return Workflow(a, b..., c..., d, f, g, h..., i..., j, l)
    end
end

function Base.read(filename::AbstractString, ::Type{Recipe})
    dict = load(filename)
    type = pop!(dict, "recipe")
    return Recipe(type, dict)
end

end
