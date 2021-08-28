module QuasiHarmonicApproxWorkflow

using AbInitioSoftwareBase: save, load
using ..Express: Calculation

export QuasiHarmonicApprox, MakeInput, CalculateThermodyn, Plot, buildworkflow

struct QuasiHarmonicApprox <: Calculation end

include("Config.jl")

function buildjob end

function buildworkflow(cfgfile)
    return chain(
        buildjob(MakeInput{QuasiHarmonicApprox}(), cfgfile),
        buildjob(CalculateThermodyn{QuasiHarmonicApprox}(), cfgfile),
        buildjob(Plot{QuasiHarmonicApprox}(), cfgfile),
    )
end

include("DefaultActions.jl")
using .DefaultActions: MakeInput, CalculateThermodyn, Plot

include("Recipes.jl")

end
