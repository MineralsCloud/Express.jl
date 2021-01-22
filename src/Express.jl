module Express

using AbInitioSoftwareBase: load
using AbInitioSoftwareBase.Inputs: Input
using AbInitioSoftwareBase.CLI: Mpiexec, scriptify
using Mustache: render_from_file, render, @mt_str
using SimpleWorkflow: Script, ExternalAtomicJob, parallel
using Unitful: uparse
import Unitful
import UnitfulAtomic

myuparse(string) = uparse(string; unit_context = [Unitful, UnitfulAtomic])

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
    return quotient,
    Tuple(range(quotient * i; stop = quotient * (i + 1) - 1) for i in 0:(njobs-1))
end

function makescript_from_file(to, file, view)
    str = render_from_file(file, view)
    return Script(str, to)
end
function makescript(to, template, view)
    str = render(template, view)
    return Script(str, to)
end

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

struct MakeCmd{T} <: Action{T} end
MakeCmd(::T) where {T<:Calculation} = MakeCmd{T}()
function (::MakeCmd)(
    output,
    input,
    np,
    exe;
    use_shell = false,
    script_template = nothing,
    shell_args = Dict(),
    procs = (),
    kwargs...,
)
    if isnothing(script_template)
        return scriptify(
            Mpiexec(np, Pair[]),
            exe;
            stdin = input,
            stdout = output,
            use_shell = use_shell,
        )
    else
        view = merge(
            Dict(
                "output" => output,
                "input" => input,
                "np" => np,
                "exe" => exe,
                "script_template" => script_template,
                "procs" => procs,
            ),
            shell_args,
        )
        saveto, _ = mktemp(dirname(input); cleanup = false)
        return makescript_from_file(saveto, script_template, view)
    end
end
function (x::MakeCmd)(outputs::AbstractArray, inputs::AbstractArray, np, exe; kwargs...)
    # `map` guarantees they are of the same size, no need to check.
    n, proc_sets = distprocs(np, length(inputs))
    return map(outputs, inputs, proc_sets) do output, input, procs
        x(output, input, n, exe; procs = procs, kwargs...)
    end
end

function buildjob(x::MakeCmd, outputs, inputs, np, exe; kwargs...)
    jobs = map(ExternalAtomicJob, x(outputs, inputs, np, exe; kwargs...))
    return parallel(jobs...)
end

function buildworkflow(file)
    config = load(file)
    mod = whichmodule(config["workflow"])
    return mod.buildworkflow(file)
end

function loadconfig(file)
    config = load(file)
    push!(config, "workdir" => abspath(dirname(file)))  # Add `workdir` key since we now deprecate it
    mod = whichmodule(config["workflow"])
    mod.checkconfig(config)  # Errors will be thrown if exist
    return mod.materialize(config)
end

function currentsoftware end

# include("SelfConsistentField.jl")
# include("BandStructure.jl")
include("EosFitting/EosFitting.jl")
include("Phonon/Phonon.jl")
include("Qha/Qha.jl")

end
