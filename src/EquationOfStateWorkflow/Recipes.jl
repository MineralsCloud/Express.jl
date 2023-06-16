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

function stage(::Scf, r::ParallelEosFittingRecipe)
    steps = map((
        DownloadPotentials{Scf}(),
        MakeInput{Scf}(),
        RunCmd{Scf}(),
        ExtractData{Scf}(),
        SaveVolumeEnergy{Scf}(),
        FitEquationOfState{Scf}(),
        SaveParameters{Scf}(),
    )) do action
        think(action, r.config)
    end
    download = Job(steps[1]; name="download potentials")
    makeinputs = map(thunk -> Job(thunk; name="make input in SCF"), steps[2])
    runcmds = map(
        thunk -> WeaklyDependentJob(thunk; name="run ab initio software in SCF"), steps[3]
    )
    extractdata = map(
        thunk -> WeaklyDependentJob(thunk; name="extract E(V) data in SCF"), steps[4]
    )
    savedata = StronglyDependentJob(steps[5]; name="save E(V) data in SCF")
    fiteos = StronglyDependentJob(steps[6]; name="fit E(V) data in SCF")
    saveparams = StronglyDependentJob(steps[7]; name="save EOS parameters in SCF")
    download .→ makeinputs .→ runcmds .→ extractdata .→ fiteos → saveparams
    extractdata .→ savedata
    return steps = (;
        download=download,
        makeinputs=makeinputs,
        runcmds=runcmds,
        extractdata=extractdata,
        savedata=savedata,
        fiteos=fiteos,
        saveparams=saveparams,
    )
end
function stage(::VariableCellOptimization, r::ParallelEosFittingRecipe)
    steps = map((
        MakeInput{VariableCellOptimization}(),
        RunCmd{VariableCellOptimization}(),
        ExtractData{VariableCellOptimization}(),
        SaveVolumeEnergy{VariableCellOptimization}(),
        FitEquationOfState{VariableCellOptimization}(),
        SaveParameters{VariableCellOptimization}(),
    )) do action
        think(action, r.config)
    end
    makeinputs = map(
        thunk -> StronglyDependentJob(thunk; name="make input in vc-relax"), steps[1]
    )
    runcmds = map(
        thunk -> WeaklyDependentJob(thunk; name="run ab initio software in vc-relax"),
        steps[2],
    )
    extractdata = map(
        thunk -> WeaklyDependentJob(thunk; name="extract E(V) data in vc-relax"), steps[3]
    )
    savedata = StronglyDependentJob(steps[4]; name="save E(V) data in vc-relax")
    fiteos = StronglyDependentJob(steps[5]; name="fit E(V) data in vc-relax")
    saveparams = StronglyDependentJob(steps[6]; name="save EOS parameters in vc-relax")
    makeinputs .→ runcmds .→ extractdata .→ fiteos → saveparams
    extractdata .→ savedata
    return steps = (;
        makeinputs=makeinputs,
        runcmds=runcmds,
        extractdata=extractdata,
        savedata=savedata,
        fiteos=fiteos,
        saveparams=saveparams,
    )
end

function build(::Type{Workflow}, r::ParallelEosFittingRecipe)
    stage₁, stage₂ = stage(Scf(), r), stage(VariableCellOptimization(), r)
    stage₁.fiteos .→ stage₂.makeinputs
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
