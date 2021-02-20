module Express

using AbInitioSoftwareBase: load
using AbInitioSoftwareBase.Inputs: Input
using SimpleWorkflow: Script, ExternalAtomicJob, parallel
using Unitful: uparse
import Unitful
import UnitfulAtomic

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

macro action(type)
    quote
        struct $type{T} <: Action{T} end
    end |> esc
end

calculation(::Action{T}) where {T} = T()

function whichmodule(name)
    name = lowercase(name)
    return if name == "eos"
        EosFitting
    elseif name in ("phonon dispersion", "vdos")
        Phonon
    elseif name in ("qha single", "qha multi")
        Qha
    else
        error("workflow `$name` is not recognized!")
    end
end

function buildworkflow(file)
    config = load(file)
    mod = whichmodule(config["workflow"])
    return mod.buildworkflow(file)
end

function loadconfig(file)
    config = load(file)
    if haskey(config, "dirs")
        if !haskey(config, "root")
            config["dirs"]["root"] = abspath(dirname(file))
        end
    else
        config["dirs"] = Dict("root" => abspath(dirname(file)))
    end
    mod = whichmodule(pop!(config, "workflow"))
    return mod.Config.materialize(config)
end

function currentsoftware end

# include("SelfConsistentField.jl")
# include("BandStructure.jl")
include("Shell.jl")
using .Shell
include("EosFitting/EosFitting.jl")
include("Phonon/Phonon.jl")
include("Qha/Qha.jl")

end
