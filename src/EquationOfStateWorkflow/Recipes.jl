module Recipes

using Configurations: from_dict
using EasyJobsBase: →
using ExpressBase: Scf, VariableCellOptimization
using ExpressBase.Files: load
using ExpressBase.Recipes: Recipe
using SimpleWorkflows: Workflow

using ...Express: DownloadPotentials, RunCmd, jobify
using ..Config: RuntimeConfig
using ..EquationOfStateWorkflow:
    MakeInput, SaveVolumeEnergy, FitEquationOfState, SaveParameters

struct ParallelEosFittingRecipe <: Recipe
    config
end

function build(::Type{Workflow}, r::ParallelEosFittingRecipe)
    for T in (Scf, VariableCellOptimization)
        steps = map((
            DownloadPotentials{T}(),
            MakeInput{T}(),
            RunCmd{T}(),
            SaveVolumeEnergy{T}(),
            FitEquationOfState{T}(),
            SaveParameters{T}(),
        )) do action
            jobify(action, r.config)
        end
        download, makeinput, runcmd, savedata, fiteos, saveparams = steps
        download .→ makeinput .→ runcmd .→ savedata .→ fiteos → saveparams
    end
    return Workflow(download)
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
