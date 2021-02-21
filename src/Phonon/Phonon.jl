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
using SimpleWorkflow: run!, →
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
using ..Shell: distprocs, @intjob

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
    settings = load(cfgfile)
    x = if settings["workflow"] == "phonon dispersion"
        PhononDispersion
    elseif settings["workflow"] == "vdos"
        VDos
    else
        error("unsupported option!")
    end
    return begin
        @intjob(LogMsg{Scf}()(true)) →
        @intjob(MakeInput{Scf}()(cfgfile)) →
        buildjob(MakeCmd{Scf}(), cfgfile) →
        @intjob(LogMsg{Scf}()(false)) →
        @intjob(LogMsg{Dfpt}()(true)) →
        @intjob(MakeInput{Dfpt}()(cfgfile)) →
        buildjob(MakeCmd{Dfpt}(), cfgfile)[1] →
        @intjob(LogMsg{Dfpt}()(false)) →
        @intjob(LogMsg{RealSpaceForceConstants}()(true)) →
        @intjob(MakeInput{RealSpaceForceConstants}()(cfgfile)) →
        buildjob(MakeCmd{RealSpaceForceConstants}(), cfgfile) →
        @intjob(LogMsg{RealSpaceForceConstants}()(false)) →
        @intjob(LogMsg{x}()(true)) →
        @intjob(MakeInput{x}()(cfgfile)) →
        buildjob(MakeCmd{x}(), cfgfile)[1] → @intjob(LogMsg{x}()(false))
    end
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
