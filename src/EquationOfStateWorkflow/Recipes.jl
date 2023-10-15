module Recipes

using Configurations: from_dict
using EasyJobsBase: Job, ConditionalJob, ArgDependentJob, →
using ExpressBase: SelfConsistentField, VariableCellOptimization, think
using ExpressBase.Files: load
using ExpressBase.Recipes: Recipe
using SimpleWorkflows: Workflow, eachjob, run!

using ..Config: StaticConfig, expand
using ..EquationOfStateWorkflow:
    DownloadPotentials,
    ComputeVolume,
    CreateInput,
    WriteInput,
    RunCmd,
    ExtractData,
    GatherData,
    SaveData,
    FitEquationOfState,
    SaveParameters

export build, run!

struct ParallelEosFittingRecipe <: Recipe
    config
end

function stage(::SelfConsistentField, r::ParallelEosFittingRecipe)
    conf = expand(r.config, SelfConsistentField())
    steps = map((
        DownloadPotentials(SelfConsistentField()),
        ComputeVolume(SelfConsistentField()),
        CreateInput(SelfConsistentField()),
        WriteInput(SelfConsistentField()),
        RunCmd(SelfConsistentField()),
        ExtractData(SelfConsistentField()),
        GatherData(SelfConsistentField()),
        SaveData(SelfConsistentField()),
        FitEquationOfState(SelfConsistentField()),
        SaveParameters(SelfConsistentField()),
    )) do action
        think(action, conf)
    end
    download = Job(steps[1]; name="download potentials")
    compute = map(thunk -> Job(thunk; name="compute volume in SCF"), steps[2])
    makeinputs = map(thunk -> ArgDependentJob(thunk; name="update input in SCF"), steps[3])
    writeinputs = map(thunk -> ArgDependentJob(thunk; name="write input in SCF"), steps[4])
    runcmds = map(
        thunk -> ConditionalJob(thunk; name="run ab initio software in SCF"), steps[5]
    )
    extractdata = map(
        thunk -> ConditionalJob(thunk; name="extract E(V) data in SCF"), steps[6]
    )
    gatherdata = ArgDependentJob(steps[7]; name="gather E(V) data in SCF")
    savedata = ArgDependentJob(steps[8]; name="save E(V) data in SCF")
    fiteos = ArgDependentJob(steps[9]; name="fit E(V) data in SCF")
    saveparams = ArgDependentJob(steps[10]; name="save EOS parameters in SCF")
    download .→
    compute .→
    makeinputs .→ writeinputs .→ runcmds .→ extractdata .→ gatherdata → fiteos → saveparams
    gatherdata → savedata
    return steps = (;
        download=download,
        compute=compute,
        makeinputs=makeinputs,
        writeinputs=writeinputs,
        runcmds=runcmds,
        extractdata=extractdata,
        gatherdata=gatherdata,
        savedata=savedata,
        fiteos=fiteos,
        saveparams=saveparams,
    )
end
function stage(::VariableCellOptimization, r::ParallelEosFittingRecipe)
    conf = expand(r.config, VariableCellOptimization())
    steps = map((
        ComputeVolume(VariableCellOptimization()),
        CreateInput(VariableCellOptimization()),
        WriteInput(VariableCellOptimization()),
        RunCmd(VariableCellOptimization()),
        ExtractData(VariableCellOptimization()),
        GatherData(VariableCellOptimization()),
        SaveData(VariableCellOptimization()),
        FitEquationOfState(VariableCellOptimization()),
        SaveParameters(VariableCellOptimization()),
    )) do action
        think(action, conf)
    end
    compute = map(thunk -> Job(thunk; name="compute volume in vc-relax"), steps[1])
    makeinputs = map(
        thunk -> ArgDependentJob(thunk; name="make input in vc-relax"), steps[2]
    )
    writeinputs = map(
        thunk -> ArgDependentJob(thunk; name="write input in vc-relax"), steps[3]
    )
    runcmds = map(
        thunk -> ConditionalJob(thunk; name="run ab initio software in vc-relax"), steps[4]
    )
    extractdata = map(
        thunk -> ConditionalJob(thunk; name="extract E(V) data in vc-relax"), steps[5]
    )
    gatherdata = ArgDependentJob(steps[6]; name="gather E(V) data in vc-relax")
    savedata = ArgDependentJob(steps[7]; name="save E(V) data in vc-relax")
    fiteos = ArgDependentJob(steps[8]; name="fit E(V) data in vc-relax")
    saveparams = ArgDependentJob(steps[9]; name="save EOS parameters in vc-relax")
    compute .→
    makeinputs .→ writeinputs .→ runcmds .→ extractdata .→ gatherdata → fiteos → saveparams
    gatherdata → savedata
    return steps = (;
        compute=compute,
        makeinputs=makeinputs,
        writeinputs=writeinputs,
        runcmds=runcmds,
        extractdata=extractdata,
        gatherdata=gatherdata,
        savedata=savedata,
        fiteos=fiteos,
        saveparams=saveparams,
    )
end

function build(::Type{Workflow}, r::ParallelEosFittingRecipe)
    stage₁, stage₂ = stage(SelfConsistentField(), r), stage(VariableCellOptimization(), r)
    stage₁.saveparams .→ stage₂.compute
    return Workflow(stage₁.download)
end
function build(::Type{Workflow}, file)
    dict = load(file)
    recipe = dict["recipe"]
    config = from_dict(StaticConfig, dict)
    if recipe == "eos"
        return build(Workflow, ParallelEosFittingRecipe(config))
    else
        error("unsupported recipe $recipe.")
    end
end

end
