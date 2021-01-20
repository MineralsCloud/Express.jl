module Qha

using AbInitioSoftwareBase: save, load
using SimpleWorkflow: chain

using ..Express: Calculation

export QuasiHarmonicApprox, MakeInput, CalculateThermodyn, Plot, buildworkflow

struct QuasiHarmonicApprox <: Calculation end

include("config.jl")

function buildjob end

function buildworkflow(cfgfile)
    return chain(
        buildjob(MakeInput{QuasiHarmonicApprox}(), cfgfile),
        buildjob(CalculateThermodyn{QuasiHarmonicApprox}(), cfgfile),
        buildjob(Plot{QuasiHarmonicApprox}(), cfgfile),
    )
end

module DefaultActions

using AbInitioSoftwareBase.Inputs: Input, writetxt
using PyQHA: converter, runcode, plot
using SimpleWorkflow: InternalAtomicJob

using ...Express: Action, loadconfig
using ..Qha: QuasiHarmonicApprox, materialize
import ..Qha: buildjob

include("MakeInput.jl")
include("Calculate.jl")
include("Plot.jl")

end

using .DefaultActions: MakeInput, CalculateThermodyn, Plot

end
