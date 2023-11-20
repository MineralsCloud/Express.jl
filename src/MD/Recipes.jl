module Recipes

using Configurations: from_dict
using EasyJobsBase: Job, ConditionalJob, ArgDependentJob, →
using ExpressBase:
    FixedCellOptimization,
    VariableCellOptimization,
    IonDynamics,
    VariableCellMolecularDynamics,
    think
using ExpressBase.Files: load
using ExpressBase.Recipes: Recipe
using SimpleWorkflows: Workflow, run!

using ..Config: StaticConfig, expand
using ..MD: DownloadPotentials, CreateInput, WriteInput, RunCmd, ExtractCell

struct IonDynamicsRecipe <: Recipe
    config
end
struct VariableCellMolecularDynamicsRecipe <: Recipe
    config
end

function stage(::FixedCellOptimization, r::IonDynamicsRecipe)
    conf = expand(r.config, FixedCellOptimization())
    steps = map((
        DownloadPotentials(FixedCellOptimization()),
        ExtractCell(FixedCellOptimization()),
    )) do action
        think(action, conf)
    end
    steps = Iterators.Stateful(steps)
    download = Job(first(iterate(steps)); name="download potentials")
    extractcells = map(
        thunk -> Job(thunk; name="extract cell in relax"), first(iterate(steps))
    )
    download .→ extractcells
    return (; download, extractcells)
end
function stage(::VariableCellOptimization, r::VariableCellMolecularDynamicsRecipe)
    conf = expand(r.config, VariableCellOptimization())
    steps = map((
        DownloadPotentials(VariableCellOptimization()),
        ExtractCell(VariableCellOptimization()),
    )) do action
        think(action, conf)
    end
    steps = Iterators.Stateful(steps)
    download = Job(first(iterate(steps)); name="download potentials")
    extractcells = map(
        thunk -> Job(thunk; name="extract cell in vc-relax"), first(iterate(steps))
    )
    download .→ extractcells
    return (; download, extractcells)
end
function stage(::IonDynamics, r::IonDynamicsRecipe)
    conf = expand(r.config, IonDynamics())
    steps = map((
        DownloadPotentials(IonDynamics()),
        CreateInput(IonDynamics()),
        WriteInput(IonDynamics()),
        RunCmd(IonDynamics()),
        # ExtractData(IonDynamics()),
        # GatherData(IonDynamics()),
        # SaveData(IonDynamics()),
    )) do action
        think(action, conf)
    end
    steps = Iterators.Stateful(steps)
    makeinputs = map(
        thunk -> ArgDependentJob(thunk; name="update input in MD"), first(iterate(steps))
    )
    writeinputs = map(
        thunk -> ArgDependentJob(thunk; name="write input in MD"), first(iterate(steps))
    )
    runcmds = map(
        thunk -> ConditionalJob(thunk; name="run ab initio software in MD"),
        first(iterate(steps)),
    )
    # extractdata = map(
    #     thunk -> ConditionalJob(thunk; name="extract E(V) data in MD"),
    #     first(iterate(steps)),
    # )
    # gatherdata = ArgDependentJob(first(iterate(steps)); name="gather E(V) data in MD")
    # savedata = ArgDependentJob(first(iterate(steps)); name="save E(V) data in MD")
    makeinputs .→ writeinputs .→ runcmds
    return steps = (;
        makeinputs,
        writeinputs,
        runcmds,
        # extractdata,
        # gatherdata,
        # savedata,
    )
end
function stage(::VariableCellMolecularDynamics, r::VariableCellMolecularDynamicsRecipe)
    conf = expand(r.config, VariableCellMolecularDynamics())
    steps = map((
        CreateInput(VariableCellMolecularDynamics()),
        WriteInput(VariableCellMolecularDynamics()),
        RunCmd(VariableCellMolecularDynamics()),
        # ExtractData(VariableCellMolecularDynamics()),
        # GatherData(VariableCellMolecularDynamics()),
        # SaveData(VariableCellMolecularDynamics()),
    )) do action
        think(action, conf)
    end
    steps = Iterators.Stateful(steps)
    makeinputs = map(
        thunk -> ArgDependentJob(thunk; name="update input in MD"), first(iterate(steps))
    )
    writeinputs = map(
        thunk -> ArgDependentJob(thunk; name="write input in MD"), first(iterate(steps))
    )
    runcmds = map(
        thunk -> ConditionalJob(thunk; name="run ab initio software in MD"),
        first(iterate(steps)),
    )
    # extractdata = map(
    #     thunk -> ConditionalJob(thunk; name="extract E(V) data in MD"),
    #     first(iterate(steps)),
    # )
    # gatherdata = ArgDependentJob(first(iterate(steps)); name="gather E(V) data in MD")
    # savedata = ArgDependentJob(first(iterate(steps)); name="save E(V) data in MD")
    makeinputs .→ writeinputs .→ runcmds
    return steps = (;
        makeinputs,
        writeinputs,
        runcmds,
        # extractdata,
        # gatherdata,
        # savedata,
    )
end

function build(::Type{Workflow}, r::IonDynamicsRecipe)
    stage₁ = stage(FixedCellOptimization(), r)
    stage₂ = stage(IonDynamicsRecipe(), r)
    stage₁.extractcells .→ stage₂.makeinputs
    return Workflow(stage₂.download)
end
function build(::Type{Workflow}, r::VariableCellMolecularDynamicsRecipe)
    stage₁ = stage(VariableCellOptimization(), r)
    stage₂ = stage(VariableCellMolecularDynamics(), r)
    stage₁.extractcells .→ stage₂.makeinputs
    return Workflow(stage₂.download)
end
function build(::Type{Workflow}, file)
    dict = load(file)
    recipe = dict["recipe"]
    config = from_dict(StaticConfig, dict)
    if recipe == "md"
        return build(Workflow, IonDynamicsRecipe(config))
    elseif recipe == "vc-md"
        return build(Workflow, VariableCellMolecularDynamicsRecipe(config))
    else
        error("unsupported recipe $recipe.")
    end
end

end
