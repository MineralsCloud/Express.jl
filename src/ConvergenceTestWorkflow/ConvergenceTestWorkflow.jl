module ConvergenceTestWorkflow

using ..Express: Scf, calculation

function isconvergent(a::AbstractVector)
    terms = abs.(diff(a))
    x, y, z = last(terms, 3)
    return all(0 <= r < 1 for r in (y / x, z / y))
end

include("Config.jl")
include("actions.jl")
include("Recipes.jl")

end
