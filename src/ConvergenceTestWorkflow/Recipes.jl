module Recipes

using AbInitioSoftwareBase: load
using Serialization: deserialize
using SimpleWorkflows: AtomicJob, Workflow, run!, ▷, ⋲, ⋺

using ..ConvergenceTestWorkflow:
    Scf, DownloadPotentials, LogMsg, MakeInput, RunCmd, TestConvergence, buildjob

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
        d = buildjob(TestConvergence{Scf}(), cfgfile)
        e = AtomicJob(() -> LogMsg{Scf}()(; start = false))
        ((((a0 ▷ a) ⋲ b) ▷ c) ⋺ d) ▷ e
        return Workflow(a0, a, b..., c..., d, e)
    end
end

end
