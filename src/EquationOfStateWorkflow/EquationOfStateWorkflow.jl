module EquationOfStateWorkflow

import ..Express
using ..Express:
    Optimization,
    SelfConsistentField,
    Scf,
    FixedIonSelfConsistentField,
    calculation,
    current_software

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
    RunCmd,
    LogMsg,
    calculation

struct StructuralOptimization <: Optimization end
struct VariableCellOptimization <: Optimization end
# See https://www.quantum-espresso.org/Doc/pw_user_guide/node10.html
const FixedCellOptimization = StructuralOptimization
const StOptim = StructuralOptimization
const VcOptim = VariableCellOptimization
const ScfOrOptim = Union{SelfConsistentField,Optimization}

include("Config.jl")
include("actions.jl")
include("Recipes.jl")

end
