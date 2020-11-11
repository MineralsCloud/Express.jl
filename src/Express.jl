module Express

using AbInitioSoftwareBase: load
using AbInitioSoftwareBase.CLI: Mpiexec
using Mustache: render
using SimpleWorkflow: Script, ExternalAtomicJob, InternalAtomicJob, chain, parallel

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

function whichmodule(name)
    name = lowercase(name)
    return if name == "eos"
        EosFitting
    elseif name in ("phonon dispersion", "vdos")
        Phonon
    end
end

function buildjob(outputs, inputs, np, exe; kwargs...)
    # `map` guarantees they are of the same size, no need to check.
    n = distprocs(np, length(inputs))
    subjobs = map(outputs, inputs) do output, input
        f = Mpiexec(n; kwargs...) ∘ exe
        cmd = f(stdin = input, stdout = output)
        ExternalAtomicJob(cmd)
    end
    return parallel(subjobs...)
end
function buildjob(cfgfile)
    settings = load(cfgfile)
    mod = whichmodule(settings["workflow"])
    return getproperty(mod, :buildjob)(cfgfile)
end

function buildworkflow(cfgfile)
    settings = load(cfgfile)
    mod = whichmodule(settings["workflow"])
    return getproperty(mod, :buildworkflow)(cfgfile)
end

function load_settings(cfgfile)
    settings = load(cfgfile)
    mod = whichmodule(settings["workflow"])
    getproperty(mod, :check_settings)(settings)  # Errors will be thrown if exist
    return getproperty(mod, :expand_settings)(settings)
end

# include("SelfConsistentField.jl")
# include("BandStructure.jl")
include("EosFitting.jl")
include("Phonon.jl")

end
