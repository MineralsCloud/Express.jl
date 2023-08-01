module ExpressBase

export SelfConsistentField,
    NonSelfConsistentField,
    BandStructure,
    DFTPlusU,
    FixedCellOptimization,
    VariableCellOptimization,
    DensityFunctionalPerturbationTheory,
    RealSpaceForceConstants,
    PhononDispersion,
    PhononDensityOfStates,
    ZoneCenterPhonons,
    QuasiHarmonicApproximation,
    SCF,
    NSCF,
    FixedIonSelfConsistentField,
    DFPT,
    VDOS,
    ZoneCentrePhonons,
    QHA

# See https://www.quantum-espresso.org/Doc/pw_user_guide/node10.html
"Represent all materials calculations."
abstract type Calculation end
abstract type ElectronicStructure <: Calculation end
struct SelfConsistentField <: ElectronicStructure end
struct NonSelfConsistentField <: ElectronicStructure end
struct BandStructure <: ElectronicStructure end
struct DFTPlusU <: ElectronicStructure end
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
const SCF = SelfConsistentField
const NSCF = NonSelfConsistentField
const FixedIonSelfConsistentField = SelfConsistentField
const DFPT = DensityFunctionalPerturbationTheory
const VDOS = PhononDensityOfStates
const ZoneCentrePhonons = ZoneCenterPhonons  # For British users
const QHA = QuasiHarmonicApproximation

"Represent an atomic action for a specific `Calculation` type."
abstract type Action{T<:Calculation} end

include("procs_per_job.jl")
include("Files.jl")
include("actions.jl")
include("Config.jl")
include("Recipes.jl")

end
