module Recipes

using AbInitioSoftwareBase: load
using ExpressBase: Scf
using ExpressBase.Recipes: Recipe
using ExpressWorkflowMaker.Templates: DownloadPotentials, LogTime, RunCmd, jobify
using SimpleWorkflows: Workflow, run!, →, ⇉, ⇶, ⭃

using ..ConvergenceTestWorkflow: MakeInput, TestConvergence

export build, run!

struct TestCutoffEnergyRecipe <: Recipe
    config
end
struct TestKPointsRecipe <: Recipe
    config
end

function build(::Type{Workflow}, r::TestCutoffEnergyRecipe)
    a = jobify(DownloadPotentials{Scf}(), r.config)
    b = jobify(LogTime{Scf}())
    c = jobify(MakeInput{Scf}(), r.config)
    d = jobify(RunCmd{Scf}(), r.config)
    e = jobify(TestConvergence{Scf}(), r.config)
    a → b ⇉ c ⇶ d ⭃ e
    return Workflow(a)
end
function build(::Type{Workflow}, r::TestKPointsRecipe)
    a = jobify(DownloadPotentials{Scf}(), r.config)
    b = jobify(LogTime{Scf}())
    c = jobify(MakeInput{Scf}(), r.config)
    d = jobify(RunCmd{Scf}(), r.config)
    e = jobify(TestConvergence{Scf}(), r.config)
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
