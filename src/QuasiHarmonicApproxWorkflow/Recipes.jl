module Recipes

using SimpleWorkflows: AtomicJob, Workflow, run!, ▷
using ..QuasiHarmonicApproxWorkflow: QuasiHarmonicApprox
using ..QuasiHarmonicApproxWorkflow.DefaultActions:
    MakeInput, CalculateThermodyn, Plot, buildjob

function buildworkflow(cfgfile)
    a = buildjob(MakeInput{QuasiHarmonicApprox}(), cfgfile)
    b = buildjob(CalculateThermodyn{QuasiHarmonicApprox}(), cfgfile)
    c = buildjob(Plot{QuasiHarmonicApprox}(), cfgfile)
    a ▷ b ▷ c
    return Workflow(a, b, c)
end

end
