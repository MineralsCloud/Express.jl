module Qha

using AbInitioSoftwareBase: save, load

using ..Express: Calculation

export QuasiHarmonicApprox, MakeInput, CalculateThermodyn, Plot

struct QuasiHarmonicApprox <: Calculation end

include("config.jl")

module DefaultActions

using AbInitioSoftwareBase.Inputs: Input, writeinput
using PyQHA: converter, runcode, plot
using SimpleWorkflow: InternalAtomicJob

using ...Express: Action, loadconfig
using ..Qha: QuasiHarmonicApprox, materialize

include("MakeInput.jl")
include("Calculate.jl")
include("Plot.jl")

end

using .DefaultActions: MakeInput, CalculateThermodyn, Plot

end
