module ExpressBase

export SelfConsistentField,
    FixedIonSelfConsistentField,
    NonSelfConsistentField,
    BandStructure,
    DFTPlusU,
    IonDynamics,
    VariableCellMolecularDynamics,
    FixedCellOptimization,
    VariableCellOptimization,
    LinearResponse,
    FourierTransform,
    PhononDispersion,
    PhononDensityOfStates,
    ZoneCenterPhonons,
    ZoneCentrePhonons,
    QuasiHarmonicApproximation,
    SCF,
    NSCF,
    VDOS,
    QHA

# See https://www.quantum-espresso.org/Doc/pw_user_guide/node10.html
"Represent all materials calculations."
abstract type Calculation end
abstract type ElectronicStructure <: Calculation end
struct SelfConsistentField <: ElectronicStructure end
struct NonSelfConsistentField <: ElectronicStructure end
struct BandStructure <: ElectronicStructure end
struct DFTPlusU <: ElectronicStructure end
# MD
abstract type MolecularDynamics <: Calculation end
struct IonDynamics <: MolecularDynamics end
struct VariableCellMolecularDynamics <: MolecularDynamics end
# Optimization
abstract type Optimization <: Calculation end
struct FixedCellOptimization <: Optimization end
struct VariableCellOptimization <: Optimization end
# Phonon
abstract type LatticeDynamics <: Calculation end
struct LinearResponse <: LatticeDynamics end
struct FourierTransform <: LatticeDynamics end
struct PhononDispersion <: LatticeDynamics end
struct PhononDensityOfStates <: LatticeDynamics end
struct ZoneCenterPhonons <: LatticeDynamics end
struct QuasiHarmonicApproximation <: Calculation end
# Aliases
const SCF = SelfConsistentField
const NSCF = NonSelfConsistentField
const FixedIonSelfConsistentField = SelfConsistentField
const VDOS = PhononDensityOfStates
const ZoneCentrePhonons = ZoneCenterPhonons  # For British users
const QHA = QuasiHarmonicApproximation

"Represent an atomic action for a specific `Calculation` type."
abstract type Action{T<:Calculation} end

Calculation(::Action{T}) where {T} = T()

include("procs_per_job.jl")
include("Files.jl")
include("actions.jl")
include("Config.jl")
include("Recipes.jl")

end
