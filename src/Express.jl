module Express

export Step, PWscfCalculation, PHononCalculation, CPCalculation

struct Step{N} end
Step(n) = Step{n}()

abstract type QuantumESPRESSOCalculation end
struct PWscfCalculation{T} <: QuantumESPRESSOCalculation end
struct PHononCalculation{T} <: QuantumESPRESSOCalculation end
struct CPCalculation{T} <: QuantumESPRESSOCalculation end

include("Jobs.jl")
include("SelfConsistentField.jl")
include("BandStructure.jl")
include("EquationOfStateFitting.jl")
# include("Phonon.jl")
include("Wizard/Wizard.jl")

end # module
