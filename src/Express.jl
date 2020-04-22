module Express

export Step,
    ScfCalculation,
    PhononCalculation,
    StructureOptimization,
    CPMD,
    PrepareInput,
    LaunchJob,
    AnalyseOutput

struct Step{N} end
Step(N::Integer) = N > 0 ? Step{N}() : throw(ArgumentError("step `$N` is nonpositive!"))

abstract type Calculation end
struct ScfCalculation <: Calculation end
struct StructureOptimization <: Calculation end
struct PhononCalculation <: Calculation end
struct CPMD <: Calculation end

abstract type Procedure end
struct PrepareInput <: Procedure end
struct LaunchJob <: Procedure end
struct AnalyseOutput <: Procedure end

include("CLI.jl")
include("Jobs.jl")
include("SelfConsistentField.jl")
include("BandStructure.jl")
include("EquationOfStateFitting.jl")
# include("Phonon.jl")
# include("Wizard/Wizard.jl")

end # module
