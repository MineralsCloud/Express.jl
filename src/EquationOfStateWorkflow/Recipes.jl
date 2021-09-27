module Recipes

using AbInitioSoftwareBase: load
using Serialization: deserialize
using SimpleWorkflows: AtomicJob, Workflow, run!, ▷, ⋲, ⋺

using ..Config: ExpandConfig
using ..EquationOfStateWorkflow:
    Scf, VcOptim, LogMsg, MakeInput, RunCmd, GetData, FitEos, buildjob

export buildworkflow, run!

struct Recipe{T}
    config::Dict{String,Any}
end

function buildworkflow(cfgfile)
    dict = load(cfgfile)
    if isfile(dict["save"]["status"])
        w = deserialize(dict["save"]["status"])
        typeassert(w, Workflow)
        return w
    else
        a = AtomicJob(() -> LogMsg{Scf}()(; start = true))
        b = buildjob(MakeInput{Scf}(), cfgfile)
        c = buildjob(RunCmd{Scf}(), cfgfile)
        d = buildjob(FitEos{Scf}(), cfgfile)
        f = AtomicJob(() -> LogMsg{Scf}()(; start = false))
        g = AtomicJob(() -> LogMsg{VcOptim}()(; start = true))
        h = buildjob(MakeInput{VcOptim}(), cfgfile)
        i = buildjob(RunCmd{VcOptim}(), cfgfile)
        j = buildjob(FitEos{VcOptim}(), cfgfile)
        l = AtomicJob(() -> LogMsg{VcOptim}()(; start = false))
        (((((((a ⋲ b) ▷ c) ⋺ d) ▷ f) ▷ g ⋲ h) ▷ i) ⋺ j) ▷ l
        return Workflow(a, b..., c..., d, f, g, h..., i..., j, l)
    end
end

end
