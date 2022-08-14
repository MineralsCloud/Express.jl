module ExpressBase

export SelfConsistentField,
    NonSelfConsistentField,
    BandStructure,
    DftPlusU,
    FixedCellOptimization,
    VariableCellOptimization,
    DensityFunctionalPerturbationTheory,
    RealSpaceForceConstants,
    PhononDispersion,
    PhononDensityOfStates,
    ZoneCenterPhonons,
    QuasiHarmonicApproximation,
    Scf,
    NScf,
    FixedIonSelfConsistentField,
    Dfpt,
    VDos,
    ZoneCentrePhonons,
    Qha
export calculation

# See https://www.quantum-espresso.org/Doc/pw_user_guide/node10.html
"Represent all materials calculations."
abstract type Calculation end
abstract type ElectronicStructure <: Calculation end
struct SelfConsistentField <: ElectronicStructure end
struct NonSelfConsistentField <: ElectronicStructure end
struct BandStructure <: ElectronicStructure end
struct DftPlusU <: ElectronicStructure end
# Optimization
abstract type Optimization <: Calculation end
struct FixedCellOptimization <: Optimization end
struct VariableCellOptimization <: Optimization end
# Phonon
abstract type LatticeDynamics <: Calculation end
struct DensityFunctionalPerturbationTheory <: LatticeDynamics end
struct RealSpaceForceConstants <: LatticeDynamics end
struct PhononDispersion <: LatticeDynamics end
struct PhononDensityOfStates <: LatticeDynamics end
struct ZoneCenterPhonons <: LatticeDynamics end
struct QuasiHarmonicApproximation <: Calculation end
# Aliases
const Scf = SelfConsistentField
const NScf = NonSelfConsistentField
const FixedIonSelfConsistentField = SelfConsistentField
const Dfpt = DensityFunctionalPerturbationTheory
const VDos = PhononDensityOfStates
const ZoneCentrePhonons = ZoneCenterPhonons  # For British users
const Qha = QuasiHarmonicApproximation

"Represent an atomic action for a specific `Calculation` type."
abstract type Action{T<:Calculation} end

"""
    calculation(::Action)

Return the calculation type of the `Action`.
"""
calculation(::Action{T}) where {T} = T()

include("Recipes.jl")

end
