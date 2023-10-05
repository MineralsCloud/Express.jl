module Recipes

using Configurations: from_dict
using EasyJobsBase: Job, ConditionalJob, ArgDependentJob, →
using ExpressBase: SelfConsistentField, think
using ExpressBase.Files: load
using ExpressBase.Recipes: Recipe
using SimpleWorkflows: Workflow, eachjob, run!

using ..Config: StaticConfig, expand
using ..ConvergenceTestWorkflow:
    DownloadPotentials,
    CreateInput,
    WriteInput,
    RunCmd,
    ExtractData,
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
        SaveData(SelfConsistentField()),
        TestConvergence(SelfConsistentField()),
    )) do action
        think(action, conf)
    end
    download = Job(steps[1]; name="download potentials")
    makeinputs = map(thunk -> Job(thunk; name="update input in SCF"), steps[2])
    writeinputs = map(thunk -> ArgDependentJob(thunk; name="write input in SCF"), steps[3])
    runcmds = map(
        thunk -> ConditionalJob(thunk; name="run ab initio software in SCF"), steps[4]
    )
    extractdata = map(
        thunk -> ConditionalJob(thunk; name="extract energies in SCF"), steps[5]
    )
    savedata = ArgDependentJob(steps[6]; name="save energies in SCF")
    testconv = ArgDependentJob(steps[9]; name="test convergence in SCF")
    download .→ makeinputs .→ writeinputs .→ runcmds .→ extractdata .→ testconv
    extractdata .→ savedata
    return steps = (;
        download=download,
        makeinputs=makeinputs,
        writeinputs=writeinputs,
        runcmds=runcmds,
        extractdata=extractdata,
        savedata=savedata,
        testconv=testconv,
    )
end

function build(::Type{Workflow}, r::Union{TestCutoffEnergyRecipe,TestKMeshRecipe})
    stage = stage(SelfConsistentField(), r)
    return Workflow(stage.download)
end
function build(::Type{Workflow}, file)
    dict = load(file)
    recipe = dict["recipe"]
    if recipe == "ecut"
        return build(Workflow, TestCutoffEnergyRecipe(dict))
    elseif recipe == "k_mesh"
        return build(Workflow, TestKMeshRecipe(dict))
    else
        error("unsupported recipe $recipe.")
    end
end

end
