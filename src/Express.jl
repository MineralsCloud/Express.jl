module Express

abstract type Calculation end
abstract type ElectronicStructure <: Calculation end
abstract type Optimization <: Calculation end
abstract type VibrationalProperty <: Calculation end

# include("SelfConsistentField.jl")
# include("BandStructure.jl")
include("EosFitting/EosFitting.jl")
include("Phonon/Phonon.jl")

end # module
