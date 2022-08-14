module Recipes

using AbInitioSoftwareBase: load
using SimpleWorkflows: Job, Workflow, run!, →
using ExpressBase: QuasiHarmonicApproximation
using ExpressBase.Recipes: Recipe
using ExpressWorkflowMaker.Templates: jobify

using ..QuasiHarmonicApproxWorkflow: MakeInput, CalculateThermodyn, Plot

export build, run!

struct SingleConfigurationRecipe <: Recipe
    config
end
struct MultiConfigurationRecipe <: Recipe
    config
end

function build(
    ::Type{Workflow},
    r::Union{SingleConfigurationRecipe,MultiConfigurationRecipe},
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
