module EosFitting

using Serialization: deserialize
using SimpleWorkflow: Workflow, run!, →

import ..Express
using ..Express:
    Optimization,
    SelfConsistentField,
    Scf,
    FixedIonSelfConsistentField,
    calculation,
    currentsoftware,
    loadconfig
using ..Shell: @intjob

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
    config = loadconfig(cfgfile)
    if isfile(config.recover)
        w = deserialize(config.recover)
        typeassert(w, Workflow)
        return w
    else
        return begin
            @intjob(LogMsg{Scf}()(true)) →
            @intjob(MakeInput{Scf}()(cfgfile)) →
            buildjob(MakeCmd{Scf}(), cfgfile) →
            @intjob(FitEos{Scf}()(cfgfile)) →
            @intjob(
                GetData{Scf}()(shortname(Scf) * ".json", last.(iofiles(Scf(), cfgfile)))
            ) →
            @intjob(LogMsg{Scf}()(false)) →
            @intjob(LogMsg{VcOptim}()(true)) →
            @intjob(MakeInput{VcOptim}()(cfgfile)) →
            buildjob(MakeCmd{VcOptim}(), cfgfile) →
            @intjob(FitEos{VcOptim}()(cfgfile)) →
            @intjob(
                GetData{VcOptim}()(
                    shortname(VcOptim) * ".json",
                    last.(iofiles(VcOptim(), cfgfile)),
                )
            ) → @intjob(LogMsg{VcOptim}()(false))
        end
    end
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

buildjob(x::Union{MakeInput,GetData,FitEos,LogMsg}, args...) =
    InternalAtomicJob(() -> x(args...))

end

using .DefaultActions: MakeInput, GetData, FitEos, SaveEos, MakeCmd, LogMsg

end
