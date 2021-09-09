module Express

using AbInitioSoftwareBase: load
using AbInitioSoftwareBase.Inputs: Input
using Unitful: uparse
import Unitful
import UnitfulAtomic

import Configurations: convert_to_option

myuparse(str::AbstractString) =
    uparse(filter(!isspace, str); unit_context = [Unitful, UnitfulAtomic])
myuparse(num::Number) = num  # FIXME: this might be error-prone!

abstract type Calculation end
abstract type ElectronicStructure <: Calculation end
struct SelfConsistentField <: ElectronicStructure end
abstract type Optimization <: Calculation end
abstract type LatticeDynamics <: Calculation end
# Aliases
const Calc = Calculation
const Optim = Optimization
const Scf = SelfConsistentField
const FixedIonSelfConsistentField = SelfConsistentField

abstract type Action{T<:Calculation} end

calculation(::Action{T}) where {T} = T()

function whichmodule(name)
    name = lowercase(name)
    return if name == "eos"
        EquationOfStateWorkflow
    elseif name in ("phonon dispersion", "vdos")
        # PhononWorkflow
    elseif name in ("qha single", "qha multi")
        # QuasiHarmonicApproxWorkflow
    else
        error("workflow `$name` is not recognized!")
    end
end

function buildworkflow(file)
    config = load(file)
    mod = whichmodule(config["workflow"])
    return mod.buildworkflow(file)
end
abstract type UnitfulVector end

function current_software end
convert_to_option(::Type{<:UnitfulVector}, ::Type{AbstractVector}, s) = eval(Meta.parse(s))

# include("SelfConsistentField.jl")
# include("BandStructure.jl")
include("Shell.jl")
using .Shell
include("EquationOfStateWorkflow/EquationOfStateWorkflow.jl")
include("PhononWorkflow/PhononWorkflow.jl")
include("QuasiHarmonicApproxWorkflow/QuasiHarmonicApproxWorkflow.jl")

end
