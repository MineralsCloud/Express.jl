module Express

using Mustache: render
using SimpleWorkflow: Script

abstract type Calculation end
abstract type ElectronicStructure <: Calculation end
struct SelfConsistentField <: ElectronicStructure end
abstract type Optimization <: Calculation end
abstract type VibrationalProperty <: Calculation end
# Aliases
const Calc = Calculation
const Optim = Optimization
const Scf = SelfConsistentField
const FixedIonSelfConsistentField = SelfConsistentField

function distprocs(nprocs, njobs)
    quotient, remainder = divrem(nprocs, njobs)
    if !iszero(remainder)
        @warn "The processes are not fully balanced! Consider the number of subjobs!"
    end
    return quotient
end

function makescript(template, view)
    map((:press, :nprocs, :in, :out, :script)) do key
        @assert haskey(view, key)
    end
    str = render(template, view)
    return Script(str, view[:script])
end
makescript(template, args::Pair...) = makescript(template, Dict(args))
makescript(template; kwargs...) = makescript(template, Dict(kwargs))

# include("SelfConsistentField.jl")
# include("BandStructure.jl")
include("EosFitting.jl")
include("Phonon.jl")

end # module
