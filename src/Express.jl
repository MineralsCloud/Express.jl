module Express

using Dates: DateTime, format
using Unitful
using UnitfulAtomic

abstract type Calculation end
abstract type ElectronicStructure <: Calculation end
struct SelfConsistentField <: ElectronicStructure end
struct NonSelfConsistentField <: ElectronicStructure end
struct BandStructure <: ElectronicStructure end
struct ElectronicDensityOfStates <: ElectronicStructure end
struct DftPlusU <: ElectronicStructure end
struct HartreeFock <: ElectronicStructure end
abstract type Optimization <: Calculation end
struct StructuralOptimization <: Optimization end
struct VariableCellOptimization <: Optimization end
abstract type Dynamics <: Calculation end
struct MolecularDynamics <: Dynamics end
struct VariableCellMolecularDynamics <: Dynamics end
abstract type VibrationalProperty <: Calculation end
struct DfptMethod <: VibrationalProperty end
struct SmallDisplacementMethod <: VibrationalProperty end
struct ForceConstant <: VibrationalProperty end
struct PhononDispersion <: VibrationalProperty end
struct PhononDensityOfStates <: VibrationalProperty end

struct Action{T} end

struct Step{S<:Calculation,T<:Action} end
Step(c::Calculation, a::Action) = Step{typeof(c),typeof(a)}()
(c::Calculation)(a::Action) = Step(c, a)

calctype(::Type{Step{S,T}}) where {S,T} = S
calctype(s::Step) = calctype(typeof(s))  # No instance, `S` could be abstract

actiontype(::Type{Step{S,T}}) where {S,T} = T
actiontype(s::Step) = actiontype(typeof(s))

abstract type StepStatus end
struct Pending <: StepStatus end
abstract type Finished <: StepStatus end
struct Succeeded <: Finished end
struct Failed <: Finished end

mutable struct Context
    inputs
    outputs
    status::StepStatus
    time::DateTime
    step::Step
end

_emoji(::Pending) = '🚧'
_emoji(::Succeeded) = '✅'
_emoji(::Failed) = '❌'

_uparse(str::AbstractString) = uparse(str; unit_context = [Unitful, UnitfulAtomic])

Base.show(io::IO, ::SelfConsistentField) = print(io, "scf calculation")
Base.show(io::IO, ::NonSelfConsistentField) = print(io, "nscf calculation")
Base.show(io::IO, ::StructuralOptimization) = print(io, "relax calculation")
Base.show(io::IO, ::VariableCellOptimization) = print(io, "vc-relax calculation")
Base.show(io::IO, ::MolecularDynamics) = print(io, "md calculation")
Base.show(io::IO, ::VariableCellMolecularDynamics) = print(io, "vc-md calculation")
Base.show(io::IO, ::Action{:prepare_input}) = print(io, "prepare input")
Base.show(io::IO, ::Action{:launch_job}) = print(io, "launch job")
Base.show(io::IO, ::Action{:analysis_output}) = print(io, "analysis output")
Base.show(io::IO, ::Step{S,T}) where {S,T} = print(io, "$(T()) @ $(S())")
function Base.show(io::IO, x::Context)
    print(io, _emoji(x.status), ' ')
    printstyled(io, x.step; bold = true)
    printstyled(io, " @ ", format(x.time, "Y/mm/dd H:M:S"); color = :light_black)
end # function Base.show

include("Jobs.jl")
# include("SelfConsistentField.jl")
# include("BandStructure.jl")
include("EosFitting/EosFitting.jl")
include("Phonon/Phonon.jl")
# include("Wizard/Wizard.jl")

end # module
