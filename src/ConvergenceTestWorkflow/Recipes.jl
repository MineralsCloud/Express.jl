module Recipes

using ExpressBase: SCF
using ExpressBase.Files: load
using ExpressBase.Recipes: Recipe
using SimpleWorkflows.Workflows: Workflow, run!, →, ⇉, ⇶, ⭃

using ...Express: DownloadPotentials, RunCmd, jobify
using ..ConvergenceTestWorkflow: MakeInput, TestConvergence

export build, run!

struct TestCutoffEnergyRecipe <: Recipe
    config
end
struct TestKPointsRecipe <: Recipe
    config
end

function build(::Type{Workflow}, r::TestCutoffEnergyRecipe)
    a = jobify(DownloadPotentials{SCF}(), r.config)
    c = jobify(MakeInput{SCF}(), r.config)
    d = jobify(RunCmd{SCF}(), r.config)
    e = jobify(TestConvergence{SCF}(), r.config)
    a → b ⇉ c ⇶ d ⭃ e
    return Workflow(a)
end
function build(::Type{Workflow}, r::TestKPointsRecipe)
    a = jobify(DownloadPotentials{SCF}(), r.config)
    c = jobify(MakeInput{SCF}(), r.config)
    d = jobify(RunCmd{SCF}(), r.config)
    e = jobify(TestConvergence{SCF}(), r.config)
    a → b ⇉ c ⇶ d ⭃ e
    return Workflow(a)
end
function build(::Type{Workflow}, file)
    dict = load(file)
    recipe = dict["recipe"]
    if recipe == "ecut"
        return build(Workflow, TestCutoffEnergyRecipe(dict))
    elseif recipe == "k_mesh"
        return build(Workflow, TestKPointsRecipe(dict))
    else
        error("unsupported recipe $recipe.")
    end
end

end
