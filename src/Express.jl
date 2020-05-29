module Express

using JSON
using YAML
using Unitful
using UnitfulAtomic

@enum Action begin
    PREPARE_INPUT = 1
    LAUNCH_JOB = 2
    ANALYSE_OUTPUT = 3
end

abstract type Calculation{T} end
struct SelfConsistentField{T} <: Calculation{T} end

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

include("CLI.jl")
include("Jobs.jl")
include("Schemes.jl")
# include("SelfConsistentField.jl")
# include("BandStructure.jl")
include("EosFitting.jl")
include("Phonon.jl")
# include("Wizard/Wizard.jl")

end # module
