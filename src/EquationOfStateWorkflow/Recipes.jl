module Recipes

using AbInitioSoftwareBase: load
using ExpressBase: Scf, VariableCellOptimization
using ExpressWorkflowMaker.Templates: DownloadPotentials, LogTime, RunCmd, jobify
using SimpleWorkflows: Job, Workflow, run!, →, ⇉, ⇶, ⭃

using ..Config: ExpandConfig
using ..EquationOfStateWorkflow: MakeInput, GetData, FitEos, buildjob

export buildworkflow, run!

function buildworkflow(cfgfile)
    stage = Scf
    a0 = buildjob(DownloadPotentials{stage}(), cfgfile)
    a = Job(() -> LogTime{stage}()())
    b = buildjob(MakeInput{stage}(), cfgfile)
    c = buildjob(RunCmd{stage}(), cfgfile)
    d0 = buildjob(GetData{stage}(), cfgfile)
    d = buildjob(FitEos{stage}(), cfgfile)
    stage = VariableCellOptimization
    g = Job(() -> LogTime{stage}()())
    h = buildjob(MakeInput{stage}(), cfgfile)
    i = buildjob(RunCmd{stage}(), cfgfile)
    j0 = buildjob(GetData{stage}(), cfgfile)
    j = buildjob(FitEos{stage}(), cfgfile)
    a0 → a ⇉ b ⇶ c ⭃ d0 → d → g ⇉ h ⇶ i ⭃ j0 → j
    return Workflow(a0)
end

end
