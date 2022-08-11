module EquationOfStateWorkflow

using ..Express: current_software

export MakeInput, FitEos, GetData, SaveEos, RunCmd, LogMsg, calculation

const ScfOrOptim = Union{SelfConsistentField,Optimization}

include("Config.jl")
include("actions.jl")
include("Recipes.jl")

end
