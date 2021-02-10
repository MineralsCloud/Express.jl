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
using ..Shell: MakeCmd, distprocs

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
    makescript,
    run!,
    loadconfig

struct Dfpt <: LatticeDynamics end
struct RealSpaceForceConstants <: LatticeDynamics end
struct PhononDispersion <: LatticeDynamics end
struct PhononDensityOfStates <: LatticeDynamics end
struct ZoneCenterPhonons <: LatticeDynamics end
const VDos = PhononDensityOfStates
const ZoneCentrePhonons = ZoneCenterPhonons  # For English users

function buildjob end

function buildworkflow(cfgfile)
    step1 = buildjob(MakeInput{Scf}(), cfgfile)
    step2 = chain(step1, buildjob(MakeCmd{Scf}(), cfgfile)[1])
    step3 = chain(step2[end], buildjob(MakeInput{Dfpt}(), cfgfile))
    step4 = chain(step3[end], buildjob(MakeCmd{Dfpt}(), cfgfile)[1])
    step5 = chain(step4[end], buildjob(MakeInput{RealSpaceForceConstants}(), cfgfile))
    step6 = chain(step5[end], buildjob(MakeCmd{RealSpaceForceConstants}(), cfgfile)[1])
    settings = load(cfgfile)
    x = if settings["workflow"] == "phonon dispersion"
        PhononDispersion
    elseif settings["workflow"] == "vdos"
        VDos
    else
        error("unsupported option!")
    end
    step7 = chain(step6[end], buildjob(MakeInput{x}(), cfgfile))
    step8 = chain(step7[end], buildjob(MakeCmd{x}(), cfgfile)[1])
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
using SimpleWorkflow: InternalAtomicJob

using ...Express: Action, Calculation, LatticeDynamics, Scf, loadconfig
using ...EosFitting: VcOptim
using ..Phonon:
    Dfpt, RealSpaceForceConstants, PhononDispersion, VDos, shortname, prevcalc, order

import ..Phonon: buildjob

include("MakeInput.jl")

end

using .DefaultActions: MakeInput

end
