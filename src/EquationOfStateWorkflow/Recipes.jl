module Recipes

using Configurations: from_dict
using EasyJobsBase: Job, ConditionalJob, ArgDependentJob, →
using ExpressBase: Scf, VariableCellOptimization
using ExpressBase.Files: load
using ExpressBase.Recipes: Recipe
using SimpleWorkflows: Workflow, eachjob, run!

using ...Express: DownloadPotentials, RunCmd, think
using ..Config: RuntimeConfig
using ..EquationOfStateWorkflow:
    MakeInput, ExtractData, SaveData, FitEquationOfState, SaveParameters

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
        SaveData{Scf}(),
        FitEquationOfState{Scf}(),
        SaveParameters{Scf}(),
    )) do action
        think(action, r.config)
    end
    download = Job(steps[1]; name="download potentials")
    makeinputs = map(thunk -> Job(thunk; name="make input in SCF"), steps[2])
    runcmds = map(
        thunk -> ConditionalJob(thunk; name="run ab initio software in SCF"), steps[3]
    )
    extractdata = map(
        thunk -> ConditionalJob(thunk; name="extract E(V) data in SCF"), steps[4]
    )
    savedata = ArgDependentJob(steps[5]; name="save E(V) data in SCF")
    fiteos = ArgDependentJob(steps[6]; name="fit E(V) data in SCF")
    saveparams = ArgDependentJob(steps[7]; name="save EOS parameters in SCF")
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
        SaveData{VariableCellOptimization}(),
        FitEquationOfState{VariableCellOptimization}(),
        SaveParameters{VariableCellOptimization}(),
    )) do action
        think(action, r.config)
    end
    makeinputs = map(
        thunk -> ArgDependentJob(thunk; name="make input in vc-relax"), steps[1]
    )
    runcmds = map(
        thunk -> ConditionalJob(thunk; name="run ab initio software in vc-relax"), steps[2]
    )
    extractdata = map(
        thunk -> ConditionalJob(thunk; name="extract E(V) data in vc-relax"), steps[3]
    )
    savedata = ArgDependentJob(steps[4]; name="save E(V) data in vc-relax")
    fiteos = ArgDependentJob(steps[5]; name="fit E(V) data in vc-relax")
    saveparams = ArgDependentJob(steps[6]; name="save EOS parameters in vc-relax")
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
    return Workflow(stage₁.download)
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
