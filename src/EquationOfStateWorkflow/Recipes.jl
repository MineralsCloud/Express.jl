module Recipes

using Configurations: from_dict
using EasyJobs: →, ⇉, ⇶, ⭃
using ExpressBase: Scf, VariableCellOptimization
using ExpressBase.Files: load
using ExpressBase.Recipes: Recipe
using SimpleWorkflows: Workflow, run!

using ...Express: DownloadPotentials, LogTime, RunCmd, jobify
using ..Config: RuntimeConfig
using ..EquationOfStateWorkflow: MakeInput, GetRawData, FitEos

struct ParallelEosFittingRecipe <: Recipe
    config
end

function build(::Type{Workflow}, r::ParallelEosFittingRecipe)
    a = jobify(DownloadPotentials{Scf}(), r.config)
    b = jobify(LogTime{Scf}())
    c = jobify(MakeInput{Scf}(), r.config)
    d = jobify(RunCmd{Scf}(), r.config)
    e = jobify(GetRawData{Scf}(), r.config)
    f = jobify(FitEos{Scf}(), r.config)
    g = jobify(LogTime{VariableCellOptimization}())
    h = jobify(MakeInput{VariableCellOptimization}(), r.config)
    i = jobify(RunCmd{VariableCellOptimization}(), r.config)
    j = jobify(GetRawData{VariableCellOptimization}(), r.config)
    k = jobify(FitEos{VariableCellOptimization}(), r.config)
    a → b ⇉ c ⇶ d ⭃ e → f → g ⇉ h ⇶ i ⭃ j → k
    return Workflow(a)
end
function build(::Type{Workflow}, file)
    dict = load(file)
    recipe = dict["recipe"]
    config = from_dict(RuntimeConfig, dict)
    if recipe == "eos"
        return build(Workflow, ParallelEosFittingRecipe(config))
    else
        error("unsupported recipe $recipe.")
    end
end

end
