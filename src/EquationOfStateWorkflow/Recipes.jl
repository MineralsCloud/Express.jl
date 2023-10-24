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
    ExtractCell,
    SaveCell,
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
    steps = Iterators.Stateful((
        DownloadPotentials(SelfConsistentField()),
        ComputeVolume(SelfConsistentField()),
        CreateInput(SelfConsistentField()),
        WriteInput(SelfConsistentField()),
        RunCmd(SelfConsistentField()),
        ExtractCell(SelfConsistentField()),
        SaveCell(SelfConsistentField()),
        ExtractData(SelfConsistentField()),
        GatherData(SelfConsistentField()),
        SaveData(SelfConsistentField()),
        FitEquationOfState(SelfConsistentField()),
        SaveParameters(SelfConsistentField()),
    )) do action
        think(action, conf)
    end
    download = Job(iterate(steps); name="download potentials")
    compute = map(thunk -> Job(thunk; name="compute volume in SCF"), iterate(steps))
    makeinputs = map(
        thunk -> ArgDependentJob(thunk; name="update input in SCF"), iterate(steps)
    )
    writeinputs = map(
        thunk -> ArgDependentJob(thunk; name="write input in SCF"), iterate(steps)
    )
    runcmds = map(
        thunk -> ConditionalJob(thunk; name="run ab initio software in SCF"), iterate(steps)
    )
    extractcells = map(
        thunk -> ConditionalJob(thunk; name="extract cell in SCF"), iterate(steps)
    )
    savecells = map(thunk -> ConditionalJob(thunk; name="save cell in SCF"), iterate(steps))
    extractdata = map(
        thunk -> ConditionalJob(thunk; name="extract E(V) data in SCF"), iterate(steps)
    )
    gatherdata = ArgDependentJob(iterate(steps); name="gather E(V) data in SCF")
    savedata = ArgDependentJob(iterate(steps); name="save E(V) data in SCF")
    fiteos = ArgDependentJob(iterate(steps); name="fit E(V) data in SCF")
    saveparams = ArgDependentJob(iterate(steps); name="save EOS parameters in SCF")
    download .→
    compute .→
    makeinputs .→ writeinputs .→ runcmds .→ extractdata .→ gatherdata → fiteos → saveparams
    gatherdata → savedata
    runcmds .→ extractcells .→ savecells
    return steps = (;
        download=download,
        compute=compute,
        makeinputs=makeinputs,
        writeinputs=writeinputs,
        runcmds=runcmds,
        extractcells=extractcells,
        savecells=savecells,
        extractdata=extractdata,
        gatherdata=gatherdata,
        savedata=savedata,
        fiteos=fiteos,
        saveparams=saveparams,
    )
end
function stage(::VariableCellOptimization, r::ParallelEosFittingRecipe)
    conf = expand(r.config, VariableCellOptimization())
    steps = Iterators.Stateful((
        ComputeVolume(VariableCellOptimization()),
        CreateInput(VariableCellOptimization()),
        WriteInput(VariableCellOptimization()),
        RunCmd(VariableCellOptimization()),
        ExtractCell(VariableCellOptimization()),
        SaveCell(VariableCellOptimization()),
        ExtractData(VariableCellOptimization()),
        GatherData(VariableCellOptimization()),
        SaveData(VariableCellOptimization()),
        FitEquationOfState(VariableCellOptimization()),
        SaveParameters(VariableCellOptimization()),
    )) do action
        think(action, conf)
    end
    compute = map(thunk -> Job(thunk; name="compute volume in vc-relax"), iterate(steps))
    makeinputs = map(
        thunk -> ArgDependentJob(thunk; name="make input in vc-relax"), iterate(steps)
    )
    writeinputs = map(
        thunk -> ArgDependentJob(thunk; name="write input in vc-relax"), iterate(steps)
    )
    runcmds = map(
        thunk -> ConditionalJob(thunk; name="run ab initio software in vc-relax"),
        iterate(steps),
    )
    extractcells = map(
        thunk -> ConditionalJob(thunk; name="extract cell in SCF"), iterate(steps)
    )
    savecells = map(thunk -> ConditionalJob(thunk; name="save cell in SCF"), iterate(steps))
    extractdata = map(
        thunk -> ConditionalJob(thunk; name="extract E(V) data in vc-relax"), iterate(steps)
    )
    gatherdata = ArgDependentJob(iterate(steps); name="gather E(V) data in vc-relax")
    savedata = ArgDependentJob(iterate(steps); name="save E(V) data in vc-relax")
    fiteos = ArgDependentJob(iterate(steps); name="fit E(V) data in vc-relax")
    saveparams = ArgDependentJob(iterate(steps); name="save EOS parameters in vc-relax")
    compute .→
    makeinputs .→ writeinputs .→ runcmds .→ extractdata .→ gatherdata → fiteos → saveparams
    gatherdata → savedata
    runcmds .→ extractcells .→ savecells
    return steps = (;
        compute=compute,
        makeinputs=makeinputs,
        writeinputs=writeinputs,
        runcmds=runcmds,
        extractcells=extractcells,
        savecells=savecells,
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
