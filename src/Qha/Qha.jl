module Qha

using AbInitioSoftwareBase: save, load

using ..Express: Calculation

export QuasiHarmonicApprox, MakeInput, CalculateThermodyn

struct QuasiHarmonicApprox <: Calculation end

include("config.jl")

module DefaultActions

using AbInitioSoftwareBase.Inputs: Input, writeinput
using PyQHA: converter, runcode
using SimpleWorkflow: InternalAtomicJob

using ...Express: Action, loadconfig
using ..Qha: QuasiHarmonicApprox, materialize

include("MakeInput.jl")
include("Calculate.jl")

end

using .DefaultActions: MakeInput, CalculateThermodyn

end
