module Recipes

using SimpleWorkflows: AtomicJob, Workflow, run!, ▷

function buildworkflow(cfgfile)
    return buildjob(MakeInput{QuasiHarmonicApprox}(), cfgfile) ▷
           buildjob(CalculateThermodyn{QuasiHarmonicApprox}(), cfgfile) ▷
           buildjob(Plot{QuasiHarmonicApprox}(), cfgfile)
end

end
