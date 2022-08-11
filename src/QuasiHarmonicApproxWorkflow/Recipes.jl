module Recipes

using SimpleWorkflows: Job, Workflow, run!, →
using ExpressBase: QuasiHarmonicApproximation

using ..QuasiHarmonicApproxWorkflow: MakeInput, CalculateThermodyn, Plot, buildjob

export buildworkflow, run!

function buildworkflow(cfgfile)
    a = buildjob(MakeInput{QuasiHarmonicApproximation}(), cfgfile)
    b = buildjob(CalculateThermodyn{QuasiHarmonicApproximation}(), cfgfile)
    c = buildjob(Plot{QuasiHarmonicApproximation}(), cfgfile)
    a → b → c
    return Workflow(a)
end

end
