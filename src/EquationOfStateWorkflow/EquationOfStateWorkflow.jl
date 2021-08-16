module EquationOfStateWorkflow

import ..Express
using ..Express:
    Optimization,
    SelfConsistentField,
    Scf,
    FixedIonSelfConsistentField,
    calculation,
    current_software,
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
    calculation

struct StructuralOptimization <: Optimization end
struct VariableCellOptimization <: Optimization end
# See https://www.quantum-espresso.org/Doc/pw_user_guide/node10.html
const FixedCellOptimization = StructuralOptimization
const StOptim = StructuralOptimization
const VcOptim = VariableCellOptimization
const ScfOrOptim = Union{SelfConsistentField,Optimization}

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
using Unitful: ustrip, unit

using ...Express: Action, loadconfig
using ..EquationOfStateWorkflow: ScfOrOptim, Scf

struct MakeCmd{T} <: Action{T} end

include("MakeInput.jl")
include("GetData.jl")
include("FitEos.jl")
include("SaveEos.jl")
include("LogMsg.jl")

end

using .DefaultActions: MakeInput, GetData, FitEos, SaveEos, MakeCmd, LogMsg

include("Recipes.jl")

end
