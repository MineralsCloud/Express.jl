module Recipes

using AbInitioSoftwareBase: load
using ExpressBase: Scf
using ExpressWorkflowMaker.Templates: DownloadPotentials, LogTime, RunCmd, jobify
using SimpleWorkflows: Job, Workflow, run!, →, ⇉, ⇶, ⭃

using ..ConvergenceTestWorkflow: MakeInput, TestConvergence

export buildworkflow, run!

function buildworkflow(cfgfile)
    a0 = jobify(DownloadPotentials{Scf}(), cfgfile)
    a = jobify(LogTime{Scf}())
    b = jobify(MakeInput{Scf}(), cfgfile)
    c = jobify(RunCmd{Scf}(), cfgfile)
    d = jobify(TestConvergence{Scf}(), cfgfile)
    a0 → a ⇉ b ⇶ c ⭃ d
    return Workflow(a0)
end

end
