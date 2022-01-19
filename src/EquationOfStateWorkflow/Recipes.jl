module Recipes

using AbInitioSoftwareBase: load
using Serialization: deserialize
using SimpleWorkflows: AtomicJob, Workflow, run!, ▷, ⋲, ⋺

using ..Config: ExpandConfig
using ..EquationOfStateWorkflow:
    Scf, VcOptim, DownloadPotentials, LogMsg, MakeInput, RunCmd, GetData, FitEos, buildjob

export buildworkflow, run!

function buildworkflow(cfgfile)
    dict = load(cfgfile)
    if isfile(dict["save"]["status"])
        w = deserialize(dict["save"]["status"])
        typeassert(w, Workflow)
        return w
    else
        a0 = buildjob(DownloadPotentials{Scf}(), cfgfile)
        a = AtomicJob(() -> LogMsg{Scf}()(; start = true))
        b = buildjob(MakeInput{Scf}(), cfgfile)
        c = buildjob(RunCmd{Scf}(), cfgfile)
        d0 = buildjob(GetData{Scf}(), cfgfile)
        d = buildjob(FitEos{Scf}(), cfgfile)
        f = AtomicJob(() -> LogMsg{Scf}()(; start = false))
        g = AtomicJob(() -> LogMsg{VcOptim}()(; start = true))
        h = buildjob(MakeInput{VcOptim}(), cfgfile)
        i = buildjob(RunCmd{VcOptim}(), cfgfile)
        j0 = buildjob(GetData{VcOptim}(), cfgfile)
        j = buildjob(FitEos{VcOptim}(), cfgfile)
        l = AtomicJob(() -> LogMsg{VcOptim}()(; start = false))
        ((((((((((a0 ▷ a) ⋲ b) ▷ c) ⋺ d0) ▷ d) ▷ f) ▷ g ⋲ h) ▷ i) ⋺ j0) ▷ j) ▷ l
        return Workflow(a0, a, b..., c..., d0, d, f, g, h..., i..., j0, j, l)
    end
end

end
