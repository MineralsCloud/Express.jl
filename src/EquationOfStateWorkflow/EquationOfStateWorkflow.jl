module EquationOfStateWorkflow

using ExpressWorkflowMaker.Templates: LogTime, RunCmd

using ..Express: current_software

export MakeInput, FitEos, GetData, SaveEos, RunCmd, LogTime, calculation

const ScfOrOptim = Union{SelfConsistentField,Optimization}

include("Config.jl")
include("actions.jl")
include("Recipes.jl")

end
