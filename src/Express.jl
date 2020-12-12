module Express

using AbInitioSoftwareBase: load
using AbInitioSoftwareBase.Inputs: Input
using AbInitioSoftwareBase.CLI: Mpiexec, scriptify
using Mustache: render
using SimpleWorkflow: Script, ExternalAtomicJob, parallel

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
    else
        error("workflow `$name` is not recognized!")
    end
end

struct MakeCmd{T} <: Action{T} end
MakeCmd(::T) where {T<:Calculation} = MakeCmd{T}()
function (::MakeCmd)(output, input, np, exe; use_shell = false, kwargs...)
    return scriptify(
        Mpiexec(np; kwargs...),
        exe;
        stdin = input,
        stdout = output,
        use_shell = use_shell,
    )
end
function (x::MakeCmd)(outputs::AbstractArray, inputs::AbstractArray, np, exe; kwargs...)
    # `map` guarantees they are of the same size, no need to check.
    n = distprocs(np, length(inputs))
    return map(outputs, inputs) do output, input
        x(output, input, n, exe; kwargs...)
    end
end

function buildjob(x::MakeCmd, outputs, inputs, np, exe; kwargs...)
    jobs = map(ExternalAtomicJob, x(outputs, inputs, np, exe; kwargs...))
    return parallel(jobs...)
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
include("EosFitting/EosFitting.jl")
# include("Phonon.jl")

end
