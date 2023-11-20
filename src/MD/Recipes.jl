module Recipes

using Configurations: from_dict
using EasyJobsBase: Job, ConditionalJob, ArgDependentJob, →
using ExpressBase: IonDynamics, VariableCellMolecularDynamics, think
using ExpressBase.Files: load
using ExpressBase.Recipes: Recipe
using SimpleWorkflows: Workflow, run!

using ..Config: StaticConfig, expand

struct IonDynamicsRecipe <: Recipe
    config
end
struct VariableCellMolecularDynamicsRecipe <: Recipe
    config
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
    download = Job(first(iterate(steps)); name="download potentials")
    makeinputs = map(thunk -> Job(thunk; name="update input in MD"), first(iterate(steps)))
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
    download .→ makeinputs .→ writeinputs .→ runcmds
    return steps = (;
        download=download,
        makeinputs=makeinputs,
        writeinputs=writeinputs,
        runcmds=runcmds,
        # extractdata=extractdata,
        # gatherdata=gatherdata,
        # savedata=savedata,
    )
end
function stage(::VariableCellMolecularDynamics, r::VariableCellMolecularDynamicsRecipe)
    conf = expand(r.config, VariableCellMolecularDynamics())
    steps = map((
        DownloadPotentials(VariableCellMolecularDynamics()),
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
    download = Job(first(iterate(steps)); name="download potentials")
    makeinputs = map(thunk -> Job(thunk; name="update input in MD"), first(iterate(steps)))
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
    download .→ makeinputs .→ writeinputs .→ runcmds
    return steps = (;
        download=download,
        makeinputs=makeinputs,
        writeinputs=writeinputs,
        runcmds=runcmds,
        # extractdata=extractdata,
        # gatherdata=gatherdata,
        # savedata=savedata,
    )
end

function build(::Type{Workflow}, r::IonDynamicsRecipe)
    stages = [stage(IonDynamics(), r)]
    return Workflow(stages[1].download)
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
