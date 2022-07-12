module Recipes

using AbInitioSoftwareBase: load
using SimpleWorkflows: Job, Workflow, run!, →, ⇉, ⇶, ⭃

using ..ConvergenceTestWorkflow:
    Scf, DownloadPotentials, LogMsg, MakeInput, RunCmd, TestConvergence, buildjob

export buildworkflow, run!

function buildworkflow(cfgfile)
    a0 = buildjob(DownloadPotentials{Scf}(), cfgfile)
    a = Job(() -> LogMsg{Scf}()(; start = true))
    b = buildjob(MakeInput{Scf}(), cfgfile)
    c = buildjob(RunCmd{Scf}(), cfgfile)
    d = buildjob(TestConvergence{Scf}(), cfgfile)
    e = Job(() -> LogMsg{Scf}()(; start = false))
    a0 → a ⇉ b ⇶ c ⭃ d → e
    return Workflow(a0)
end

end
