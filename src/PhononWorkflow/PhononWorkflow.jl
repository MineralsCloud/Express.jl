module PhononWorkflow

using AbInitioSoftwareBase: load
using AbInitioSoftwareBase.Inputs: Input
using Serialization: deserialize
using SimpleWorkflows: Workflow, run!, ▷
using Unitful: ustrip, @u_str

import ..Express
using ..Express:
    Calculation,
    LatticeDynamics,
    SelfConsistentField,
    Scf,
    FixedIonSelfConsistentField,
    Action,
    current_software,
    myuparse
using ..EquationOfStateWorkflow: VcOptim
using ..Shell: distprocs

export Dfpt,
    Dfpt,
    SelfConsistentField,
    Scf,
    FixedIonSelfConsistentField,
    RealSpaceForceConstants,
    PhononDispersion,
    PhononDensityOfStates,
    VDos,
    ZoneCenterPhonons,
    ZoneCentrePhonons,
    MakeInput,
    LogMsg,
    run!,
    buildworkflow,
    buildjob

struct LinearResponse <: LatticeDynamics end
struct RealSpaceForceConstants <: LatticeDynamics end
struct PhononDispersion <: LatticeDynamics end
struct PhononDensityOfStates <: LatticeDynamics end
struct ZoneCenterPhonons <: LatticeDynamics end
const Dfpt = LinearResponse
const VDos = PhononDensityOfStates
const ZoneCentrePhonons = ZoneCenterPhonons  # For English users

include("Config.jl")
include("DefaultActions.jl")
using .DefaultActions: MakeInput, RunCmd, LogMsg

include("Recipes.jl")

end
