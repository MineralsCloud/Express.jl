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
using SimpleWorkflow: chain
using Unitful: ustrip, @u_str

import ..Express
using ..Express:
    Calculation,
    LatticeDynamics,
    SelfConsistentField,
    Scf,
    FixedIonSelfConsistentField,
    Action,
    MakeCmd,
    distprocs,
    currentsoftware,
    makescript,
    loadconfig,
    myuparse
using ..EosFitting: VcOptim

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
    loadconfig

struct Dfpt <: LatticeDynamics end
struct RealSpaceForceConstants <: LatticeDynamics end
struct PhononDispersion <: LatticeDynamics end
struct PhononDensityOfStates <: LatticeDynamics end
struct ZoneCenterPhonons <: LatticeDynamics end
const VDos = PhononDensityOfStates
const ZoneCentrePhonons = ZoneCenterPhonons  # For English users

function buildjob(x::MakeCmd{T}, cfgfile) where {T}
    settings = loadconfig(cfgfile)
    inp = map(dir -> joinpath(dir, shortname(T) * ".in"), settings.dirs)
    out = map(dir -> joinpath(dir, shortname(T) * ".out"), settings.dirs)
    return Express.buildjob(
        x,
        out,
        inp,
        settings.manager.np,
        settings.bin[order(T)];
        use_shell = settings.use_shell,
    )
end

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

include("config.jl")

module DefaultActions

using AbInitioSoftwareBase.Inputs: Input, writeinput
using SimpleWorkflow: InternalAtomicJob

using ...Express: Action, Calculation, LatticeDynamics, Scf, loadconfig
using ...EosFitting: VcOptim
using ..Phonon:
    Dfpt, RealSpaceForceConstants, PhononDispersion, VDos, shortname, prevcalc, order

include("makeinput.jl")

end

using .DefaultActions: MakeInput

end
