module EosFitting

using SimpleWorkflow: ExternalAtomicJob, chain, parallel, run!

import ..Express
using ..Express:
    Optimization,
    SelfConsistentField,
    Scf,
    FixedIonSelfConsistentField,
    calculation,
    currentsoftware,
    makescript,
    loadconfig

export SelfConsistentField,
    Scf,
    FixedIonSelfConsistentField,
    StructuralOptimization,
    FixedCellOptimization,
    VariableCellOptimization,
    StOptim,
    VcOptim,
    MakeInput,
    FitEos,
    GetData,
    SaveEos,
    MakeCmd,
    LogMsg,
    loadconfig,
    iofiles,
    calculation,
    makescript,
    run!,
    buildworkflow,
    buildjob

struct StructuralOptimization <: Optimization end
struct VariableCellOptimization <: Optimization end
# See https://www.quantum-espresso.org/Doc/pw_user_guide/node10.html
const FixedCellOptimization = StructuralOptimization
const StOptim = StructuralOptimization
const VcOptim = VariableCellOptimization
const ScfOrOptim = Union{SelfConsistentField,Optimization}

function iofiles(T::ScfOrOptim, cfgfile)
    config = loadconfig(cfgfile)
    return map(config.dirs) do dir
        prefix = joinpath(dir, shortname(T))
        prefix * ".in" => prefix * ".out"
    end
end

function buildjob end

function buildworkflow(cfgfile)
    step0 = buildjob(LogMsg{Scf}(), true)
    step1 = buildjob(MakeInput{Scf}(), cfgfile)
    step1 = chain(step0, step1)
    step12 = chain(step1[end], buildjob(MakeCmd{Scf}(), cfgfile)[1])
    step123 = chain(step12[end], buildjob(FitEos{Scf}(), cfgfile))
    step3p =
        buildjob(GetData{Scf}(), shortname(Scf) * ".json", last.(iofiles(Scf(), cfgfile)))
    step3p = chain(step3p, buildjob(LogMsg{Scf}(), false))
    step123 = chain(step123[end], step3p[1])
    step4m1 = buildjob(LogMsg{VcOptim}(), true)
    step4 = buildjob(MakeInput{VcOptim}(), cfgfile)
    step4 = chain(step4m1, step4)
    step45 = chain(step4[end], buildjob(MakeCmd{VcOptim}(), cfgfile)[1])
    step456 = chain(step45[end], buildjob(FitEos{VcOptim}(), cfgfile))
    step6p = buildjob(
        GetData{VcOptim}(),
        shortname(VcOptim) * ".json",
        last.(iofiles(VcOptim(), cfgfile)),
    )
    step6p = chain(step6p, buildjob(LogMsg{VcOptim}(), false))
    step456 = chain(step456[end], step6p[1])
    step16 = chain(step123[end], step456[1])
    return step16
end

shortname(calc::ScfOrOptim) = shortname(typeof(calc))

include("Config.jl")

module DefaultActions

using AbInitioSoftwareBase: save, load, extension
using AbInitioSoftwareBase.Inputs: Input, writetxt
using Dates: now, format
using EquationsOfStateOfSolids:
    EquationOfStateOfSolids, EnergyEquation, PressureEquation, Parameters, getparam
using EquationsOfStateOfSolids.Fitting: eosfit
using Logging: with_logger, current_logger
using Serialization: serialize, deserialize
using SimpleWorkflow: InternalAtomicJob
using Unitful: ustrip, unit

using ...Express: Action, loadconfig, @action
using ..EosFitting: ScfOrOptim, Scf, iofiles, shortname
import ...EosFitting: buildjob

@action MakeCmd

include("MakeInput.jl")
include("GetData.jl")
include("FitEos.jl")
include("SaveEos.jl")
include("LogMsg.jl")

end

using .DefaultActions: MakeInput, GetData, FitEos, SaveEos, MakeCmd, LogMsg

end
