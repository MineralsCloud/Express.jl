module Recipes

using ExpressBase: QuasiHarmonicApproximation
using ExpressBase.Files: load
using ExpressBase.Recipes: Recipe
using EasyJobsBase: Job
using SimpleWorkflows: Workflow, run!, →

using ..QuasiHarmonicApproxWorkflow: MakeInput, CalculateThermodyn, Plot

import ...Express: jobify

struct SingleConfigurationRecipe <: Recipe
    config
end
struct MultiConfigurationRecipe <: Recipe
    config
end

function build(
    ::Type{Workflow}, r::Union{SingleConfigurationRecipe,MultiConfigurationRecipe}
)
    a = jobify(MakeInput{QuasiHarmonicApproximation}(), r.config)
    b = jobify(CalculateThermodyn{QuasiHarmonicApproximation}(), r.config)
    c = jobify(Plot{QuasiHarmonicApproximation}(), r.config)
    a → b → c
    return Workflow(a)
end
function build(::Type{Workflow}, file)
    dict = load(file)
    recipe = dict["recipe"]
    if recipe == "qha single"
        return build(Workflow, SingleConfigurationRecipe(dict))
    elseif recipe == "multi qha"
        return build(Workflow, MultiConfigurationRecipe(dict))
    else
        error("unsupported recipe $recipe.")
    end
end

end
