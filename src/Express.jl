module Express

abstract type Calculation end
abstract type ElectronicStructure <: Calculation end
abstract type Optimization <: Calculation end
abstract type VibrationalProperty <: Calculation end

# include("SelfConsistentField.jl")
# include("BandStructure.jl")
include("EosFitting.jl")
# include("Phonon.jl")

end # module
