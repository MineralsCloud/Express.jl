module Recipes

using AbInitioSoftwareBase: load
using ExpressBase: Scf, VariableCellOptimization
using ExpressBase.Recipes: Recipe
using SimpleWorkflows.Workflows: Workflow, run!, →, ⇉, ⇶, ⭃

using ...Express: DownloadPotentials, LogTime, RunCmd, jobify
using ..EquationOfStateWorkflow: MakeInput, GetData, FitEos

struct ParallelEosFittingRecipe <: Recipe
    config
end

function build(::Type{Workflow}, r::ParallelEosFittingRecipe)
    a = jobify(DownloadPotentials{Scf}(), r.config)
    b = jobify(LogTime{Scf}())
    c = jobify(MakeInput{Scf}(), r.config)
    d = jobify(RunCmd{Scf}(), r.config)
    e = jobify(GetData{Scf}(), r.config)
    f = jobify(FitEos{Scf}(), r.config)
    g = jobify(LogTime{VariableCellOptimization}())
    h = jobify(MakeInput{VariableCellOptimization}(), r.config)
    i = jobify(RunCmd{VariableCellOptimization}(), r.config)
    j = jobify(GetData{VariableCellOptimization}(), r.config)
    k = jobify(FitEos{VariableCellOptimization}(), r.config)
    a → b ⇉ c ⇶ d ⭃ e → f → g ⇉ h ⇶ i ⭃ j → k
    return Workflow(a)
end
function build(::Type{Workflow}, file)
    dict = load(file)
    recipe = dict["recipe"]
    if recipe == "eos"
        return build(Workflow, ParallelEosFittingRecipe(dict))
    else
        error("unsupported recipe $recipe.")
    end
end

end
