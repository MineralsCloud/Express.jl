module Recipes

using SimpleWorkflows: Job, Workflow, run!, →
using ExpressBase: QuasiHarmonicApproximation
using ExpressWorkflowMaker.Templates: jobify

using ..QuasiHarmonicApproxWorkflow: MakeInput, CalculateThermodyn, Plot

export buildworkflow, run!

function buildworkflow(cfgfile)
    a = jobify(MakeInput{QuasiHarmonicApproximation}(), cfgfile)
    b = jobify(CalculateThermodyn{QuasiHarmonicApproximation}(), cfgfile)
    c = jobify(Plot{QuasiHarmonicApproximation}(), cfgfile)
    a → b → c
    return Workflow(a)
end

end
