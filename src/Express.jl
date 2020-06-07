module Express

using Unitful
using UnitfulAtomic

abstract type Property end
struct Electronic <: Property end
struct Phononic <: Property end

abstract type Simulation end
struct SelfConsistentField <: Simulation end
struct NonSelfConsistentField <: Simulation end
struct Dispersion{T<:Property} <: Simulation end
struct DensityOfStates{T<:Property} <: Simulation end
abstract type StructureOptimization <: Simulation end
struct Relaxation <: StructureOptimization end
struct VariableCellRelaxation <: StructureOptimization end
struct MolecularDynamics <: Simulation end
struct VariableCellMolecularDynamics <: Simulation end
const BandStructure = Dispersion{Electronic}

abstract type Action end
struct PreparePotential <: Action end
struct PrepareInput <: Action end
struct LaunchJob <: Action end
struct AnalyseOutput <: Action end

struct Step{S<:Simulation,T<:Action} end

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
