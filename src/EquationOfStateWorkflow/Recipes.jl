module Recipes

using AbInitioSoftwareBase: load
using ExpressBase: Scf, VariableCellOptimization
using ExpressWorkflowMaker.Templates: DownloadPotentials, LogTime, RunCmd, jobify
using SimpleWorkflows: Job, Workflow, run!, →, ⇉, ⇶, ⭃

using ..EquationOfStateWorkflow: MakeInput, GetData, FitEos

export buildworkflow, run!

function buildworkflow(cfgfile)
    stage = Scf
    a0 = jobify(DownloadPotentials{stage}(), cfgfile)
    a = jobify(LogTime{stage}())
    b = jobify(MakeInput{stage}(), cfgfile)
    c = jobify(RunCmd{stage}(), cfgfile)
    d0 = jobify(GetData{stage}(), cfgfile)
    d = jobify(FitEos{stage}(), cfgfile)
    stage = VariableCellOptimization
    g = a = jobify(LogTime{stage}())
    h = jobify(MakeInput{stage}(), cfgfile)
    i = jobify(RunCmd{stage}(), cfgfile)
    j0 = jobify(GetData{stage}(), cfgfile)
    j = jobify(FitEos{stage}(), cfgfile)
    a0 → a ⇉ b ⇶ c ⭃ d0 → d → g ⇉ h ⇶ i ⭃ j0 → j
    return Workflow(a0)
end

end
