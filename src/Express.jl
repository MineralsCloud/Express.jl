module Express

export Step, PWscfCalculation, PHononCalculation, CPCalculation

struct Step{N} end
Step(N::Integer) = N > 0 ? Step{N}() : throw(ArgumentError("step `$N` is nonpositive!"))

abstract type QuantumESPRESSOCalculation end
struct PWscfCalculation{T} <: QuantumESPRESSOCalculation end
struct PHononCalculation{T} <: QuantumESPRESSOCalculation end
struct CPCalculation{T} <: QuantumESPRESSOCalculation end

include("CLI.jl")
include("Jobs.jl")
include("SelfConsistentField.jl")
include("BandStructure.jl")
include("EquationOfStateFitting.jl")
# include("Phonon.jl")
# include("Wizard/Wizard.jl")

end # module
