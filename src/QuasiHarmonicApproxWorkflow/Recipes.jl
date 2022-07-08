module Recipes

using SimpleWorkflows: Job, Workflow, run!, →
using ..QuasiHarmonicApproxWorkflow:
    QuasiHarmonicApprox, MakeInput, CalculateThermodyn, Plot, buildjob

export buildworkflow, run!

function buildworkflow(cfgfile)
    a = buildjob(MakeInput{QuasiHarmonicApprox}(), cfgfile)
    b = buildjob(CalculateThermodyn{QuasiHarmonicApprox}(), cfgfile)
    c = buildjob(Plot{QuasiHarmonicApprox}(), cfgfile)
    return Workflow(a, b, c)
    a → b → c
end

end
