module Express

struct Step{N} end
Step(N::Integer) = N > 0 ? Step{N}() : throw(ArgumentError("step `$N` is nonpositive!"))

abstract type Calculation end
struct ScfCalculation <: Calculation end

abstract type Action end
struct PrepareInput <: Action end
struct LaunchJob <: Action end
struct AnalyseOutput <: Action end
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

include("CLI.jl")
include("Jobs.jl")
include("Schemes.jl")
# include("SelfConsistentField.jl")
# include("BandStructure.jl")
include("EosFitting.jl")
include("Phonon.jl")
# include("Wizard/Wizard.jl")

end # module
