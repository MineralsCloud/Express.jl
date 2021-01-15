module Qha

using AbInitioSoftwareBase: save, load

using ..Express: Calculation

export MakeInput, QuasiHarmonicApprox

struct QuasiHarmonicApprox <: Calculation end

include("config.jl")

module DefaultActions

using AbInitioSoftwareBase.Inputs: Input, writeinput
using PyQHA: converter
using SimpleWorkflow: InternalAtomicJob

using ...Express: Action, loadconfig
using ..Qha: QuasiHarmonicApprox, materialize

include("MakeInput.jl")

end

using .DefaultActions: MakeInput

end
