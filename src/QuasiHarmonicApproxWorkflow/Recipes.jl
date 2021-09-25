module Recipes

using SimpleWorkflows: AtomicJob, Workflow, run!, ▷
using ..QuasiHarmonicApproxWorkflow:
    QuasiHarmonicApprox, MakeInput, CalculateThermodyn, Plot, buildjob

export buildworkflow, run!

function buildworkflow(cfgfile)
    a = buildjob(MakeInput{QuasiHarmonicApprox}(), cfgfile)
    b = buildjob(CalculateThermodyn{QuasiHarmonicApprox}(), cfgfile)
    c = buildjob(Plot{QuasiHarmonicApprox}(), cfgfile)
    a ▷ b ▷ c
    return Workflow(a, b, c)
end

end
