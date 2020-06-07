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

struct Software{T} end

function save(filepath, data)
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

function load(filepath)
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

function extension(filepath)  # From https://github.com/rofinn/FilePathsBase.jl/blob/af850a4/src/path.jl#L331-L340
    name = basename(filepath)
    tokenized = split(name, '.')
    if length(tokenized) > 1
        return lowercase(tokenized[end])
    else
        return ""
    end
end

_uparse(str::AbstractString) = uparse(str; unit_context = [Unitful, UnitfulAtomic])

function _check_settings end

function Settings end

function load_settings(path)
    settings = load(path)
    _check_settings(settings)  # Errors will be thrown if exist
    return Settings(settings)
end # function load_settings

function inputstring end

include("CLI.jl")
include("Environments.jl")
include("Jobs.jl")
# include("SelfConsistentField.jl")
# include("BandStructure.jl")
include("EosFitting/EosFitting.jl")
# include("Phonon.jl")
# include("Wizard/Wizard.jl")

end # module
