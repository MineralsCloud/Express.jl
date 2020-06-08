module Express

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

abstract type Action end
struct Prepare{T} <: Action end
struct Launch{T} <: Action end
struct Analyse{T} <: Action end
const PREPARE_POTENTIAL = Prepare{:potential}()
const PREPARE_INPUT = Prepare{:input}()
const LAUNCH_JOB = Launch{:job}()
const ANALYSE_OUTPUT = Analyse{:output}()

struct Step{S<:Calculation,T<:Action} end
(::Type{S})(::T) where {S<:Calculation,T<:Action} = Step{S,T}()
Base.show(io::IO, ::Step{S,T}) where {S,T} = print(io, "step: $T for a $S calculation")

struct Workflow{T} end

struct Software{T} end

_uparse(str::AbstractString) = uparse(str; unit_context = [Unitful, UnitfulAtomic])

function inputstring end

include("fileops.jl")

function _check_settings end

function Settings end

function load_settings(path)
    settings = load(path)
    _check_settings(settings)  # Errors will be thrown if exist
    return Settings(settings)
end # function load_settings

include("CLI.jl")
include("Environments.jl")
include("Jobs.jl")
# include("SelfConsistentField.jl")
# include("BandStructure.jl")
include("EosFitting/EosFitting.jl")
# include("Phonon.jl")
# include("Wizard/Wizard.jl")

end # module
