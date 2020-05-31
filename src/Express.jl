module Express

using JSON
using YAML
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
struct PrepareInput <: Action end
struct LaunchJob <: Action end
struct AnalyseOutput <: Action end

struct Step{S<:Simulation,T<:Action} end

function save(filepath::AbstractString, data)
    ext = extension(filepath)
    if ext ∈ (".yaml", ".yml")
        YAML.write_file(expanduser(filepath), data)
    elseif ext == ".json"
        open(expanduser(filepath), "w") do io
            JSON.print(io, data)
        end
    else
        error("unknown file extension `$ext`!")
    end
end # function save

function load(filepath::AbstractString)
    ext = extension(filepath)
    if ext ∈ (".yaml", ".yml")
        return open(expanduser(filepath), "r") do io
            YAML.load(io)
        end
    elseif ext == ".json"
        return JSON.parsefile(expanduser(filepath))
    else
        error("unknown file extension `$ext`!")
    end
end # function load

extension(filepath::AbstractString) = filepath |> splitext |> last |> lowercase

_uparse(str::AbstractString) = uparse(str; unit_context = [Unitful, UnitfulAtomic])

function _check_settings end

function Settings end

function load_settings(path::AbstractString)
    settings = load(path)
    _check_settings(settings)  # Errors will be thrown if exist
    return Settings(settings)
end # function load_settings

function inputstring end

include("CLI.jl")
include("Workspaces.jl")
include("Jobs.jl")
# include("SelfConsistentField.jl")
# include("BandStructure.jl")
include("EosFitting.jl")
# include("Phonon.jl")
# include("Wizard/Wizard.jl")

end # module
