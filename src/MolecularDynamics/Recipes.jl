module Recipes

struct IonDynamicsRecipe <: Recipe
    config
end
struct VDOSRecipe <: Recipe
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
    download .→ makeinputs .→ runcmds
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

end
