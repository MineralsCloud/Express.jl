module Recipes

using AbInitioSoftwareBase: load
using ExpressBase: Scf
using ExpressWorkflowMaker.Templates: DownloadPotentials, LogTime, RunCmd, jobify
using SimpleWorkflows: Job, Workflow, run!, →, ⇉, ⇶, ⭃

using ..ConvergenceTestWorkflow: MakeInput, TestConvergence, buildjob

export buildworkflow, run!

function buildworkflow(cfgfile)
    a0 = buildjob(DownloadPotentials{Scf}(), cfgfile)
    a = Job(() -> LogTime{Scf}()())
    b = buildjob(MakeInput{Scf}(), cfgfile)
    c = buildjob(RunCmd{Scf}(), cfgfile)
    d = buildjob(TestConvergence{Scf}(), cfgfile)
    a0 → a ⇉ b ⇶ c ⭃ d
    return Workflow(a0)
end

end
