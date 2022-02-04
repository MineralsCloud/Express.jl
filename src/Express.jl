module Express

using AbInitioSoftwareBase: load
using AbInitioSoftwareBase.Inputs: Input
using Unitful: uparse
import Unitful
import UnitfulAtomic

import Configurations: convert_to_option

myuparse(str::AbstractString) =
    uparse(filter(!isspace, str); unit_context = [Unitful, UnitfulAtomic])
myuparse(num::Number) = num  # FIXME: this might be error-prone!

abstract type Calculation end
abstract type ElectronicStructure <: Calculation end
struct SelfConsistentField <: ElectronicStructure end
abstract type Optimization <: Calculation end
abstract type LatticeDynamics <: Calculation end
# Aliases
const Calc = Calculation
const Optim = Optimization
const Scf = SelfConsistentField
const FixedIonSelfConsistentField = SelfConsistentField

abstract type Action{T<:Calculation} end

calculation(::Action{T}) where {T} = T()

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
# include("QuasiHarmonicApproxWorkflow/QuasiHarmonicApproxWorkflow.jl")

end
