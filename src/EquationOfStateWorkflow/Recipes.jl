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
    MakeInput, ExtractData, SaveVolumeEnergy, FitEquationOfState, SaveParameters

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
            ExtractData{T}(),
            SaveVolumeEnergy{T}(),
            FitEquationOfState{T}(),
            SaveParameters{T}(),
        )) do action
            think(action, r.config)
        end
        download = Job(steps[1]; name="download potentials in $T")
        makeinputs = map(thunk -> Job(thunk; name="make input in $T"), steps[2])
        runcmds = map(
            thunk -> WeaklyDependentJob(thunk; name="run ab initio software in $T"),
            steps[3],
        )
        extractdata = map(
            thunk -> StronglyDependentJob(thunk; name="extract E(V) data in $T"),
            steps[4],
        )
        savedata = StronglyDependentJob(steps[5]; name="save E(V) data in $T")
        fiteos = StronglyDependentJob(steps[6]; name="fit E(V) data in $T")
        saveparams = StronglyDependentJob(steps[7]; name="save EOS parameters in $T")
        download .→ makeinputs .→ runcmds .→ extractdata .→ fiteos → saveparams
        extractdata .→ savedata
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
