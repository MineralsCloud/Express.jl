module QuasiHarmonicApproxWorkflow

using AbInitioSoftwareBase: save, load
using ..Express: Calculation

export QuasiHarmonicApprox, MakeInput, CalculateThermodyn, Plot

struct QuasiHarmonicApprox <: Calculation end

include("Config.jl")

include("DefaultActions.jl")
using .DefaultActions: MakeInput, CalculateThermodyn, Plot

include("Recipes.jl")

end
