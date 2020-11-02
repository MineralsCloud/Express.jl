module Express

abstract type Calculation end
abstract type ElectronicStructure <: Calculation end
struct SelfConsistentField <: ElectronicStructure end
abstract type Optimization <: Calculation end
abstract type VibrationalProperty <: Calculation end
# Aliases
const Calc = Calculation
const Optim = Optimization
const Scf = SelfConsistentField

# include("SelfConsistentField.jl")
# include("BandStructure.jl")
include("EosFitting.jl")
# include("Phonon.jl")

end # module
