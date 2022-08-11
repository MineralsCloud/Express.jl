module Recipes

using AbInitioSoftwareBase: load
using SimpleWorkflows: Job, Workflow, run!, →, ⇉, ⇶, ⭃

using ..Config: ExpandConfig
using ExpressBase: Scf, VariableCellOptimization
using ..EquationOfStateWorkflow:
    DownloadPotentials, LogMsg, MakeInput, RunCmd, GetData, FitEos, buildjob

export buildworkflow, run!

function buildworkflow(cfgfile)
    stage = Scf
    a0 = buildjob(DownloadPotentials{stage}(), cfgfile)
    a = Job(() -> LogMsg{stage}()(; start = true))
    b = buildjob(MakeInput{stage}(), cfgfile)
    c = buildjob(RunCmd{stage}(), cfgfile)
    d0 = buildjob(GetData{stage}(), cfgfile)
    d = buildjob(FitEos{stage}(), cfgfile)
    f = Job(() -> LogMsg{stage}()(; start = false))
    stage = VariableCellOptimization
    g = Job(() -> LogMsg{stage}()(; start = true))
    h = buildjob(MakeInput{stage}(), cfgfile)
    i = buildjob(RunCmd{stage}(), cfgfile)
    j0 = buildjob(GetData{stage}(), cfgfile)
    j = buildjob(FitEos{stage}(), cfgfile)
    l = Job(() -> LogMsg{stage}()(; start = false))
    a0 → a ⇉ b ⇶ c ⭃ d0 → d → f → g ⇉ h ⇶ i ⭃ j0 → j → l
    return Workflow(a0)
end

end
