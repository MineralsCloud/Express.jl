module Recipes

using SimpleWorkflows: AtomicJob, Workflow, run!, ▷
using ..QuasiHarmonicApproxWorkflow: QuasiHarmonicApprox
using ..QuasiHarmonicApproxWorkflow.DefaultActions:
    MakeInput, CalculateThermodyn, Plot, buildjob

function buildworkflow(cfgfile)
    return buildjob(MakeInput{QuasiHarmonicApprox}(), cfgfile) ▷
           buildjob(CalculateThermodyn{QuasiHarmonicApprox}(), cfgfile) ▷
           buildjob(Plot{QuasiHarmonicApprox}(), cfgfile)
end

end
