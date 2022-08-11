module Express

using AbInitioSoftwareBase: load
using AbInitioSoftwareBase.Inputs: Input
using Unitful: Unitful, uparse, unitmodules

import Configurations: convert_to_option

myuparse(str::AbstractString) =
    uparse(filter(!isspace, str); unit_context = push!(unitmodules, Unitful))

function current_software end

abstract type UnitfulVector end

convert_to_option(::Type{<:UnitfulVector}, ::Type{AbstractVector}, str::AbstractString) =
    eval(Meta.parse(str))

# include("SelfConsistentField.jl")
# include("BandStructure.jl")
include("Shell.jl")
using .Shell
include("ConvergenceTestWorkflow/ConvergenceTestWorkflow.jl")
include("EquationOfStateWorkflow/EquationOfStateWorkflow.jl")
include("PhononWorkflow/PhononWorkflow.jl")
include("QuasiHarmonicApproxWorkflow/QuasiHarmonicApproxWorkflow.jl")

end
