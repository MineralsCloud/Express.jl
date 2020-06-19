module Express

using AbInitioSoftwareBase: load
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

abstract type Action end
struct Prepare{T} <: Action end
struct Launch{T} <: Action end
struct Analyse{T} <: Action end
const PREPARE_POTENTIAL = Prepare{:potential}()
const PREPARE_INPUT = Prepare{:input}()
const LAUNCH_JOB = Launch{:job}()
const ANALYSE_OUTPUT = Analyse{:output}()

struct Step{S<:Calculation,T<:Action} end
Step(c::Calculation, a::Action) = Step{typeof(c),typeof(a)}()
(c::Calculation)(a::Action) = Step(c, a)
Base.show(io::IO, ::Step{S,T}) where {S,T} = print(io, "step: $(T()) @ $(S())")

calculationtype(::Type{Step{S,T}}) where {S,T} = S
calculationtype(s::Step) = calculationtype(typeof(s))  # No instance, `S` could be abstract

actiontype(::Type{Step{S,T}}) where {S,T} = T
actiontype(s::Step) = actiontype(typeof(s))

struct Workflow{T} end

_uparse(str::AbstractString) = uparse(str; unit_context = [Unitful, UnitfulAtomic])

function _check_settings end

function Settings end

function load_settings(path)
    settings = load(path)
    _check_settings(settings)  # Errors will be thrown if exist
    return Settings(settings)
end # function load_settings

Base.show(io::IO, ::SelfConsistentField) = print(io, "scf calculation")
Base.show(io::IO, ::NonSelfConsistentField) = print(io, "nscf calculation")
Base.show(io::IO, ::StructuralOptimization) = print(io, "relax calculation")
Base.show(io::IO, ::VariableCellOptimization) = print(io, "vc-relax calculation")
Base.show(io::IO, ::MolecularDynamics) = print(io, "md calculation")
Base.show(io::IO, ::VariableCellMolecularDynamics) = print(io, "vc-md calculation")
Base.show(io::IO, ::Prepare{T}) where {T} = print(io, "prepare " * lowercase(string(T)))
Base.show(io::IO, ::Launch{T}) where {T} = print(io, "launch " * lowercase(string(T)))
Base.show(io::IO, ::Analyse{T}) where {T} = print(io, "analyse " * lowercase(string(T)))

include("CLI.jl")
include("Jobs.jl")
# include("SelfConsistentField.jl")
# include("BandStructure.jl")
include("EosFitting/EosFitting.jl")
include("Phonon/Phonon.jl")
# include("Wizard/Wizard.jl")

end # module
