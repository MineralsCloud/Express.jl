module Express

abstract type Calculation end
abstract type ElectronicStructure <: Calculation end
struct SelfConsistentField <: ElectronicStructure end
struct NonSelfConsistentField <: ElectronicStructure end
struct BandStructure <: ElectronicStructure end
struct ElectronicDensityOfStates <: ElectronicStructure end
struct DftPlusU <: ElectronicStructure end
struct HartreeFock <: ElectronicStructure end
abstract type Optimization <: Calculation end
struct StructuralOptimization <: Optimization end
struct VariableCellOptimization <: Optimization end
abstract type Dynamics <: Calculation end
struct MolecularDynamics <: Dynamics end
struct VariableCellMolecularDynamics <: Dynamics end
abstract type VibrationalProperty <: Calculation end
struct DfptMethod <: VibrationalProperty end
struct SmallDisplacementMethod <: VibrationalProperty end
struct ForceConstant <: VibrationalProperty end
struct PhononDispersion <: VibrationalProperty end
struct PhononDensityOfStates <: VibrationalProperty end

Base.show(io::IO, ::SelfConsistentField) = print(io, "scf calculation")
Base.show(io::IO, ::NonSelfConsistentField) = print(io, "nscf calculation")
Base.show(io::IO, ::StructuralOptimization) = print(io, "relax calculation")
Base.show(io::IO, ::VariableCellOptimization) = print(io, "vc-relax calculation")
Base.show(io::IO, ::MolecularDynamics) = print(io, "md calculation")
Base.show(io::IO, ::VariableCellMolecularDynamics) = print(io, "vc-md calculation")

include("Jobs.jl")
# include("SelfConsistentField.jl")
# include("BandStructure.jl")
include("EosFitting/EosFitting.jl")
# include("Phonon/Phonon.jl")
# include("Wizard/Wizard.jl")

end # module
