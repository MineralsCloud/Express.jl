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
struct PreparePotential <: Action end
struct PrepareInput <: Action end
struct LaunchJob <: Action end
struct AnalyseOutput <: Action end

struct Step{S<:Calculation,T<:Action} end

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
