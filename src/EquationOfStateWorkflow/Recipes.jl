module Recipes

using AbInitioSoftwareBase: load
using Serialization: deserialize
using SimpleWorkflows: AtomicJob, Workflow, run!, ▷, ⋲, ⋺

using ..Config: ExpandConfig
using ..EquationOfStateWorkflow:
    Scf, VcOptim, LogMsg, MakeInput, RunCmd, GetData, FitEos, buildjob

export build, run!

struct Recipe{T}
    config::Dict{String,Any}
end

function build(recipe::Recipe{:eos})
    dict = load(recipe.config)
    if isfile(dict["save"]["status"])
        w = deserialize(dict["save"]["status"])
        typeassert(w, Workflow)
        return w
    else
        a = AtomicJob(() -> LogMsg{Scf}()(; start = true))
        b = buildjob(MakeInput{Scf}(), recipe.config)
        c = buildjob(RunCmd{Scf}(), recipe.config)
        d = buildjob(FitEos{Scf}(), recipe.config)
        f = AtomicJob(() -> LogMsg{Scf}()(; start = false))
        g = AtomicJob(() -> LogMsg{VcOptim}()(; start = true))
        h = buildjob(MakeInput{VcOptim}(), recipe.config)
        i = buildjob(RunCmd{VcOptim}(), recipe.config)
        j = buildjob(FitEos{VcOptim}(), recipe.config)
        l = AtomicJob(() -> LogMsg{VcOptim}()(; start = false))
        (((((((a ⋲ b) ▷ c) ⋺ d) ▷ f) ▷ g ⋲ h) ▷ i) ⋺ j) ▷ l
        return Workflow(a, b..., c..., d, f, g, h..., i..., j, l)
    end
end

end
