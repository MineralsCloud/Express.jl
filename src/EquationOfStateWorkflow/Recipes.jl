module Recipes

using Configurations: from_dict
using EasyJobsBase: →
using ExpressBase: Scf, VariableCellOptimization
using ExpressBase.Files: load
using ExpressBase.Recipes: Recipe
using SimpleWorkflows: Workflow, eachjob

using ...Express: DownloadPotentials, RunCmd, jobify
using ..Config: RuntimeConfig
using ..EquationOfStateWorkflow:
    MakeInput, SaveVolumeEnergy, FitEquationOfState, SaveParameters

struct ParallelEosFittingRecipe <: Recipe
    config
end

function build(::Type{Workflow}, r::ParallelEosFittingRecipe)
    stage₁, stage₂ = map((Scf, VariableCellOptimization)) do T
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
        Workflow(download)
    end
    stage₁ → stage₂
    return Workflow(first(eachjob(stage₁)))
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
