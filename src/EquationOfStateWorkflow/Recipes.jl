module Recipes

using Configurations: from_dict
using EasyJobsBase: Job, WeaklyDependentJob, StronglyDependentJob, →
using ExpressBase: Scf, VariableCellOptimization
using ExpressBase.Files: load
using ExpressBase.Recipes: Recipe
using SimpleWorkflows: Workflow, eachjob, run!

using ...Express: DownloadPotentials, RunCmd, think
using ..Config: RuntimeConfig
using ..EquationOfStateWorkflow:
    MakeInput, SaveVolumeEnergy, FitEquationOfState, SaveParameters

export build, run!

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
            # SaveParameters{T}(),
        )) do action
            think(action, r.config)
        end
        download = Job(steps[1]; name="download potentials")
        makeinputs = map(x -> Job(x; name="make input"), steps[2])
        runcmds = map(x -> WeaklyDependentJob(x; name="run ab initio software"), steps[3])
        savedata = StronglyDependentJob(steps[4]; name="save E(V) data")
        fiteos = StronglyDependentJob(steps[5]; name="fit E(V) data")
        download .→ makeinputs .→ runcmds .→ savedata .→ fiteos
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
