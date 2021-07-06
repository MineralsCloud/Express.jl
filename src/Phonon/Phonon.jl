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
using Serialization: deserialize
using SimpleWorkflow: Workflow, run!, →
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
    loadconfig,
    myuparse
using ..EquationOfStateWorkflows: VcOptim
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
    run!,
    loadconfig,
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

function buildjob end

function buildworkflow(cfgfile)
    config = loadconfig(cfgfile)
    if isfile(config.recover)
        w = deserialize(config.recover)
        typeassert(w, Workflow)
        return w
    else
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
            buildjob(MakeCmd{Dfpt}(), cfgfile) →
            @intjob(LogMsg{Dfpt}()(false)) →
            @intjob(LogMsg{RealSpaceForceConstants}()(true)) →
            @intjob(MakeInput{RealSpaceForceConstants}()(cfgfile)) →
            buildjob(MakeCmd{RealSpaceForceConstants}(), cfgfile) →
            @intjob(LogMsg{RealSpaceForceConstants}()(false)) →
            @intjob(LogMsg{x}()(true)) →
            @intjob(MakeInput{x}()(cfgfile)) →
            buildjob(MakeCmd{x}(), cfgfile) → @intjob(LogMsg{x}()(false))
        end
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
using ...EquationOfStateWorkflows: VcOptim
using ..Phonon:
    Dfpt, RealSpaceForceConstants, PhononDispersion, VDos, shortname, prevcalc, order

import ..Phonon: buildjob

@action MakeCmd

include("MakeInput.jl")
include("LogMsg.jl")

end

using .DefaultActions: MakeInput, MakeCmd, LogMsg

end
