"""
# module Phonon



# Examples

```jldoctest
julia>
```
"""
module Phonon

using AbInitioSoftwareBase: load
using AbInitioSoftwareBase.Inputs: Input
using SimpleWorkflow: chain, run!
using Unitful: ustrip, @u_str

import ..Express
using ..Express:
    Calculation,
    LatticeDynamics,
    SelfConsistentField,
    Scf,
    FixedIonSelfConsistentField,
    Action,
    currentsoftware,
    makescript,
    loadconfig,
    myuparse
using ..EosFitting: VcOptim
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
    makescript,
    run!,
    loadconfig,
    buildworkflow,
    buildjob

struct Dfpt <: LatticeDynamics end
struct RealSpaceForceConstants <: LatticeDynamics end
struct PhononDispersion <: LatticeDynamics end
struct PhononDensityOfStates <: LatticeDynamics end
struct ZoneCenterPhonons <: LatticeDynamics end
const VDos = PhononDensityOfStates
const ZoneCentrePhonons = ZoneCenterPhonons  # For English users

function buildjob end

function buildworkflow(cfgfile)
    step0 = buildjob(LogMsg{Scf}(), true)
    step1 = buildjob(MakeInput{Scf}(), cfgfile)
    step1 = chain(step0, step1)
    step2 = chain(step1[end], buildjob(MakeCmd{Scf}(), cfgfile)[1])
    step2 = chain(step2[end], buildjob(LogMsg{Scf}(), false))
    step2 = chain(step2[end], buildjob(LogMsg{Dfpt}(), true))
    step3 = chain(step2[end], buildjob(MakeInput{Dfpt}(), cfgfile))
    step4 = chain(step3[end], buildjob(MakeCmd{Dfpt}(), cfgfile)[1])
    step4 = chain(step4[end], buildjob(LogMsg{Dfpt}(), false))
    step4 = chain(step4[end], buildjob(LogMsg{RealSpaceForceConstants}(), true))
    step5 = chain(step4[end], buildjob(MakeInput{RealSpaceForceConstants}(), cfgfile))
    step6 = chain(step5[end], buildjob(MakeCmd{RealSpaceForceConstants}(), cfgfile)[1])
    step6 = chain(step6[end], buildjob(LogMsg{RealSpaceForceConstants}(), false))
    settings = load(cfgfile)
    x = if settings["workflow"] == "phonon dispersion"
        PhononDispersion
    elseif settings["workflow"] == "vdos"
        VDos
    else
        error("unsupported option!")
    end
    step7 = chain(step6[end], buildjob(LogMsg{x}(), true))
    step7 = chain(step7[end], buildjob(MakeInput{x}(), cfgfile))
    step8 = chain(step7[end], buildjob(MakeCmd{x}(), cfgfile)[1])
    step8 = chain(step8[end], buildjob(LogMsg{x}(), false))
    return step8
end

order(x) = order(typeof(x))
order(::Type{Scf}) = 1
order(::Type{Dfpt}) = 2
order(::Type{RealSpaceForceConstants}) = 3
order(::Type{PhononDispersion}) = 4
order(::Type{PhononDensityOfStates}) = 4

prevcalc(x) = prevcalc(typeof(x))
prevcalc(::Type{Scf}) = VcOptim()
prevcalc(::Type{Dfpt}) = Scf()
prevcalc(::Type{RealSpaceForceConstants}) = Dfpt()
prevcalc(::Type{PhononDispersion}) = RealSpaceForceConstants()
prevcalc(::Type{PhononDensityOfStates}) = RealSpaceForceConstants()

shortname(x::Calculation) = shortname(typeof(x))

include("Config.jl")

module DefaultActions

using AbInitioSoftwareBase.Inputs: Input, writetxt
using Dates: now, format
using SimpleWorkflow: InternalAtomicJob
using Logging: with_logger, current_logger

using ...Express: Action, Calculation, LatticeDynamics, Scf, loadconfig, @action
using ...EosFitting: VcOptim
using ..Phonon:
    Dfpt, RealSpaceForceConstants, PhononDispersion, VDos, shortname, prevcalc, order

import ..Phonon: buildjob

@action MakeCmd

include("MakeInput.jl")
include("LogMsg.jl")

end

using .DefaultActions: MakeInput, MakeCmd, LogMsg

end
