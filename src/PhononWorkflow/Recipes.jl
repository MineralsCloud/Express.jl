module Recipes

using Configurations: from_dict
using EasyJobsBase: Job, ConditionalJob, ArgDependentJob, getresult, eachparent, →
using ExpressBase:
    SelfConsistentField,
    DensityFunctionalPerturbationTheory,
    RealSpaceForceConstants,
    PhononDispersion,
    PhononDensityOfStates,
    think
using ExpressBase.Files: load
using ExpressBase.Recipes: Recipe
using SimpleWorkflows: Workflow, eachjob, run!
using Thinkers: Thunk, setargs!

using ..Config: StaticConfig, expand
using ..PhononWorkflow: DownloadPotentials, GatherInput, CreateInput, WriteInput, RunCmd

struct PhononDispersionRecipe <: Recipe
    config
end
struct VDOSRecipe <: Recipe
    config
end

function stage(::SelfConsistentField, r::Recipe)
    conf = expand(r.config, SelfConsistentField())
    steps = map((
        DownloadPotentials(SelfConsistentField()),
        CreateInput(SelfConsistentField()),
        WriteInput(SelfConsistentField()),
        RunCmd(SelfConsistentField()),
        # ExtractData(SelfConsistentField()),
        # GatherData(SelfConsistentField()),
        # SaveData(SelfConsistentField()),
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
        thunk -> ConditionalJob(thunk; name="extract E(V) data in SCF"), steps[5]
    )
    gatherdata = ArgDependentJob(steps[6]; name="gather E(V) data in SCF")
    savedata = ArgDependentJob(steps[7]; name="save E(V) data in SCF")
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
function stage(::DensityFunctionalPerturbationTheory, r::Recipe)
    conf = expand(r.config, DensityFunctionalPerturbationTheory())
    steps = map((
        CreateInput(DensityFunctionalPerturbationTheory()),
        WriteInput(DensityFunctionalPerturbationTheory()),
        RunCmd(DensityFunctionalPerturbationTheory()),
        # ExtractData(DensityFunctionalPerturbationTheory()),
        # GatherData(DensityFunctionalPerturbationTheory()),
        # SaveData(DensityFunctionalPerturbationTheory()),
    )) do action
        think(action, conf)
    end
    makeinputs = map(thunk -> ArgDependentJob(thunk; name="update input in DFPT"), steps[1])
    writeinputs = map(thunk -> ArgDependentJob(thunk; name="write input in DFPT"), steps[2])
    runcmds = map(
        thunk -> ConditionalJob(thunk; name="run ab initio software in DFPT"), steps[3]
    )
    # extractdata = map(
    #     thunk -> ConditionalJob(thunk; name="extract E(V) data in DFPT"), steps[4]
    # )
    # gatherdata = ArgDependentJob(steps[5]; name="gather E(V) data in DFPT")
    # savedata = ArgDependentJob(steps[6]; name="save E(V) data in DFPT")
    makeinputs .→ writeinputs .→ runcmds
    return steps = (;
        makeinputs=makeinputs,
        writeinputs=writeinputs,
        runcmds=runcmds,
        # extractdata=extractdata,
        # gatherdata=gatherdata,
        # savedata=savedata,
    )
end
function stage(::RealSpaceForceConstants, r::Recipe)
    conf = expand(r.config, RealSpaceForceConstants())
    steps = map((
        CreateInput(RealSpaceForceConstants()),
        WriteInput(RealSpaceForceConstants()),
        RunCmd(RealSpaceForceConstants()),
        # ExtractData(RealSpaceForceConstants()),
        # GatherData(RealSpaceForceConstants()),
        # SaveData(RealSpaceForceConstants()),
    )) do action
        think(action, conf)
    end
    makeinputs = map(thunk -> ArgDependentJob(thunk; name="update input in IFC"), steps[1])
    writeinputs = map(thunk -> ArgDependentJob(thunk; name="write input in IFC"), steps[2])
    runcmds = map(
        thunk -> ConditionalJob(thunk; name="run ab initio software in IFC"), steps[3]
    )
    # extractdata = map(
    #     thunk -> ConditionalJob(thunk; name="extract E(V) data in IFC"), steps[4]
    # )
    # gatherdata = ArgDependentJob(steps[5]; name="gather E(V) data in IFC")
    # savedata = ArgDependentJob(steps[6]; name="save E(V) data in IFC")
    makeinputs .→ writeinputs .→ runcmds
    return steps = (;
        makeinputs=makeinputs,
        writeinputs=writeinputs,
        runcmds=runcmds,
        # extractdata=extractdata,
        # gatherdata=gatherdata,
        # savedata=savedata,
    )
end
function stage(::PhononDispersion, r::PhononDispersionRecipe)
    conf = expand(r.config, PhononDispersion())
    steps = map((
        WriteInput(PhononDispersion()),
        RunCmd(PhononDispersion()),
        # ExtractData(PhononDispersion()),
        # GatherData(PhononDispersion()),
        # SaveData(PhononDispersion()),
    )) do action
        think(action, conf)
    end
    makeinputs = map(
        ArgDependentJob(
            Thunk(CreateInput(PhononDispersion())(conf.template));
            skip_incomplete=false,
            name="update input in phonon dispersion",
        ) for _ in Base.OneTo(length(conf.at))
    )
    foreach(makeinputs) do job
        setargs!(
            job.core,
            first(job.core.args),
            (getresult(parent) for parent in eachparent(job))...,
        )
    end
    writeinputs = map(
        thunk -> ArgDependentJob(thunk; name="write input in phonon dispersion"), steps[1]
    )
    runcmds = map(
        thunk -> ConditionalJob(thunk; name="run ab initio software in phonon dispersion"),
        steps[2],
    )
    # extractdata = map(
    #     thunk -> ConditionalJob(thunk; name="extract E(V) data in phonon dispersion"),
    #     steps[3],
    # )
    # gatherdata = ArgDependentJob(steps[4]; name="gather E(V) data in phonon dispersion")
    # savedata = ArgDependentJob(steps[5]; name="save E(V) data in phonon dispersion")
    makeinputs .→ writeinputs .→ runcmds
    return steps = (;
        makeinputs=makeinputs,
        writeinputs=writeinputs,
        runcmds=runcmds,
        # extractdata=extractdata,
        # gatherdata=gatherdata,
        # savedata=savedata,
    )
end
function stage(::PhononDensityOfStates, r::PhononDispersionRecipe)
    conf = expand(r.config, PhononDensityOfStates())
    steps = map((
        CreateInput(PhononDensityOfStates()),
        WriteInput(PhononDensityOfStates()),
        RunCmd(PhononDensityOfStates()),
        # ExtractData(PhononDensityOfStates()),
        # GatherData(PhononDensityOfStates()),
        # SaveData(PhononDensityOfStates()),
    )) do action
        think(action, conf)
    end
    makeinputs = map(
        thunk -> ArgDependentJob(thunk; name="update input in phonon dispersion"), steps[1]
    )
    writeinputs = map(
        thunk -> ArgDependentJob(thunk; name="write input in phonon dispersion"), steps[2]
    )
    runcmds = map(
        thunk -> ConditionalJob(thunk; name="run ab initio software in phonon dispersion"),
        steps[3],
    )
    # extractdata = map(
    #     thunk -> ConditionalJob(thunk; name="extract E(V) data in phonon dispersion"),
    #     steps[4],
    # )
    # gatherdata = ArgDependentJob(steps[5]; name="gather E(V) data in phonon dispersion")
    # savedata = ArgDependentJob(steps[6]; name="save E(V) data in phonon dispersion")
    makeinputs .→ writeinputs .→ runcmds
    return steps = (;
        makeinputs=makeinputs,
        writeinputs=writeinputs,
        runcmds=runcmds,
        # extractdata=extractdata,
        # gatherdata=gatherdata,
        # savedata=savedata,
    )
end

function build(::Type{Workflow}, r::PhononDispersionRecipe)
    stages = [
        stage(SelfConsistentField(), r),
        stage(DensityFunctionalPerturbationTheory(), r),
        stage(RealSpaceForceConstants(), r),
        stage(PhononDispersion(), r),
    ]
    stage[1].makeinputs .→ stage[2].makeinputs
    stage[2].makeinputs .→ stage[3].makeinputs
    stage[2].makeinputs .→ stage[4].makeinputs
    stage[3].makeinputs .→ stage[4].makeinputs
    return Workflow(stages[1].download)
end
function build(::Type{Workflow}, r::VDOSRecipe)
    stages = [
        stage(SelfConsistentField(), r),
        stage(DensityFunctionalPerturbationTheory(), r),
        stage(RealSpaceForceConstants(), r),
        stage(PhononDensityOfStates(), r),
    ]
    stage[1].makeinputs .→ stage[2].makeinputs
    stage[2].makeinputs .→ stage[3].makeinputs
    stage[2].makeinputs .→ stage[4].makeinputs
    return Workflow(stages[1].download)
end
function build(::Type{Workflow}, file)
    dict = load(file)
    recipe = dict["recipe"]
    config = from_dict(StaticConfig, dict)
    if recipe == "phonon dispersion"
        return build(Workflow, PhononDispersionRecipe(config))
    elseif dict["recipe"] == "vdos"
        return build(Workflow, VDOSRecipe(config))
    else
        error("unsupported recipe $recipe.")
    end
end

end
