module Recipes

using Configurations: from_dict
using EasyJobsBase: Job, ConditionalJob, ArgDependentJob, →
using ExpressBase: SelfConsistentField, think
using ExpressBase.Files: load
using ExpressBase.Recipes: Recipe
using SimpleWorkflows: Workflow, run!

using ..Config: StaticConfig, expand
using ..ConvergenceTestWorkflow:
    DownloadPotentials,
    CreateInput,
    WriteInput,
    RunCmd,
    ExtractData,
    GatherData,
    SaveData,
    TestConvergence

export build, run!

struct TestCutoffEnergyRecipe <: Recipe
    config
end
struct TestKMeshRecipe <: Recipe
    config
end

function stage(::SelfConsistentField, r::TestCutoffEnergyRecipe)
    conf = expand(r.config, SelfConsistentField())
    steps = map((
        DownloadPotentials(SelfConsistentField()),
        CreateInput(SelfConsistentField()),
        WriteInput(SelfConsistentField()),
        RunCmd(SelfConsistentField()),
        ExtractData(SelfConsistentField()),
        GatherData(SelfConsistentField()),
        SaveData(SelfConsistentField()),
        TestConvergence(SelfConsistentField()),
    )) do action
        think(action, conf)
    end
    steps = Iterators.Stateful(steps)
    download = Job(first(iterate(steps)); name="download potentials")
    makeinputs = map(thunk -> Job(thunk; name="update input in SCF"), first(iterate(steps)))
    writeinputs = map(
        thunk -> ArgDependentJob(thunk; name="write input in SCF"), first(iterate(steps))
    )
    runcmds = map(
        thunk -> ConditionalJob(thunk; name="run ab initio software in SCF"),
        first(iterate(steps)),
    )
    extractdata = map(
        thunk -> ConditionalJob(thunk; name="extract energies in SCF"),
        first(iterate(steps)),
    )
    gatherdata = ArgDependentJob(first(iterate(steps)); name="gather energies in SCF")
    savedata = ArgDependentJob(first(iterate(steps)); name="save energies in SCF")
    testconv = ArgDependentJob(first(iterate(steps)); name="test convergence in SCF")
    download .→ makeinputs .→ writeinputs .→ runcmds .→ extractdata .→ gatherdata → testconv
    gatherdata → savedata
    return steps = (;
        download=download,
        makeinputs=makeinputs,
        writeinputs=writeinputs,
        runcmds=runcmds,
        extractdata=extractdata,
        gatherdata=gatherdata,
        testconv=testconv,
        savedata=savedata,
    )
end

function build(::Type{Workflow}, r::Union{TestCutoffEnergyRecipe,TestKMeshRecipe})
    s = stage(SelfConsistentField(), r)
    return Workflow(s.download)
end
function build(::Type{Workflow}, file)
    dict = load(file)
    recipe = dict["recipe"]
    config = from_dict(StaticConfig, dict)
    if recipe == "ecut"
        return build(Workflow, TestCutoffEnergyRecipe(config))
    elseif recipe == "k_mesh"
        return build(Workflow, TestKMeshRecipe(config))
    else
        error("unsupported recipe $recipe.")
    end
end

end
