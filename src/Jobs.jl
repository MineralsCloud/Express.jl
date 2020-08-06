module Jobs

using AbInitioSoftwareBase.CLI: MpiExec
using Dates: DateTime, CompoundPeriod, now, canonicalize, format
using Distributed

using ..Express: Calculation

export div_nprocs,
    launchjob,
    starttime,
    stoptime,
    timecost,
    getstatus,
    isrunning,
    issucceeded,
    isfailed,
    getstdin,
    getstdout,
    getstderr,
    getresult,
    getrunning,
    getsucceeded,
    getfailed,
    peepstdin,
    peepstdout,
    peepstderr

abstract type JobStatus end
struct Running <: JobStatus end
abstract type Finished <: JobStatus end
struct Succeeded <: Finished end
struct Failed <: Finished end

abstract type AbstractJob end

mutable struct _AtomicJobTimer
    start::DateTime
    stop::DateTime
    _AtomicJobTimer() = new()
end

mutable struct _AtomicJobRef
    ref::Future
    status::JobStatus
    _AtomicJobRef() = new()
end

abstract type AtomicJob <: AbstractJob end

struct ExternalAtomicJob <: AtomicJob
    cmd::Base.AbstractCmd
    hash::UInt
    _timer::_AtomicJobTimer
    _ref::_AtomicJobRef
    ExternalAtomicJob(cmd) =
        new(cmd, hash((now(), cmd, rand(UInt))), _AtomicJobTimer(), _AtomicJobRef())
end

struct InternalAtomicJob <: AtomicJob
    fn::Function
    hash::UInt
    _timer::_AtomicJobTimer
    _ref::_AtomicJobRef
    InternalAtomicJob(fn) =
        new(fn, hash((now(), fn, rand(UInt))), _AtomicJobTimer(), _AtomicJobRef())
end

struct SerialJobs{T<:AbstractJob} <: AbstractJob
    subjobs::Vector{T}
end

struct ParallelJobs{T<:AbstractJob} <: AbstractJob
    subjobs::Vector{T}
end

Base.:∘(a::AtomicJob, b::AtomicJob...) = SerialJobs([a, b...])
Base.:∘(a::SerialJobs, b::AtomicJob...) = SerialJobs(push!(a.subjobs, b...))
Base.:∘(a::AtomicJob, b::SerialJobs) = SerialJobs(pushfirst!(b.subjobs, a))
Base.:∘(a::SerialJobs, b::SerialJobs) = SerialJobs(append!(a.subjobs, b.subjobs))
Base.:∘(a::AtomicJob, b::ParallelJobs) = SerialJobs([a, b.subjobs])
Base.:∘(a::ParallelJobs, b::AtomicJob) = SerialJobs([a.subjobs, b])
∥(a::AtomicJob, b::AtomicJob...) = ParallelJobs([a, b...])

function launchjob(cmds, interval = 3)
    subjobs = map(cmds) do cmd
        sleep(interval)
        _launch(cmd)
    end
    return ParallelJobs(vec(subjobs))
end # function launchjob
function launchjob(tracker::ParallelJobs, interval = 3)
    subjobs = map(tracker.subjobs) do subjob
        if getstatus(subjob) isa Succeeded
            subjob
        else
            sleep(interval)
            _launch(subjob)
        end
    end
    return ParallelJobs(subjobs)
end # function launchjob
function launchjob(
    outputs,
    inputs,
    np,
    softwarecmd;
    dry_run = false,
    interval = 3,
    kwargs...,
)
    # `map` guarantees they are of the same size, no need to check.
    n = div_nprocs(np, length(inputs))
    cmds = map(inputs, outputs) do input, output  # A vector of `Cmd`s
        f = MpiExec(n; kwargs...) ∘ softwarecmd
        f(stdin = input, stdout = output)
    end
    if dry_run
        @warn "the following commands will be run:"
        return cmds
    else
        return launchjob(cmds, interval)
    end
end

function _launch(cmd::Base.AbstractCmd)
    x = ExternalAtomicJob(cmd)
    x._ref.status = Running()
    x._ref.ref = @spawn begin
        x._timer.start = now()
        ref = try
            run(cmd; wait = true)  # Must wait
        catch e
            e
        end
        x._timer.stop = now()
        if ref isa Exception  # Include all cases?
            x._ref.status = Failed()
        else
            x._ref.status = Succeeded()
        end
        ref
    end
    return x
end # function _launch
function _launch(cmd::Function)
    x = InternalAtomicJob(cmd)
    x._ref.ref = @spawn begin
        x._ref.status = Running()
        x._timer.start = now()
        ref = try
            x.fn()
        catch e
            e
        end
        x._timer.stop = now()
        if ref isa Exception  # Include all cases?
            x._ref.status = Failed()
        else
            x._ref.status = Succeeded()
        end
        ref
    end
    return x
end # function _launch
_launch(x::AtomicJob) = _launch(x.cmd)
_launch(x::ParallelJobs) = map(_launch, x.subjobs)
function _launch(x::SerialJobs) end # function _launch

Base.run(x::AtomicJob) = _launch(x)

getstatus(x::AtomicJob) = x.status
getstatus(x::ParallelJobs) = map(getstatus, x.subjobs)

starttime(x::AtomicJob) = x.starttime
starttime(x::ParallelJobs) = map(starttime, x.subjobs)

stoptime(x::AtomicJob) = isrunning(x) ? nothing : x.stoptime
stoptime(x::ParallelJobs) = map(stoptime, x.subjobs)

timecost(x::AtomicJob) = isrunning(x) ? now() - x.starttime : x.stoptime - x.starttime
timecost(x::ParallelJobs) = map(timecost, x.subjobs)

getresult(x::AtomicJob) = isrunning(x) ? nothing : fetch(x.ref)
getresult(x::ParallelJobs) = map(getresult, x.subjobs)

isrunning(x::AtomicJob) = getstatus(x) === Running()

issucceeded(x::AtomicJob) = getstatus(x) === Succeeded()

isfailed(x::AtomicJob) = getstatus(x) === Failed()

getrunning(x::ParallelJobs) = ParallelJobs(_selectby(x, Running()))

getsucceeded(x::ParallelJobs) = ParallelJobs(_selectby(x, Succeeded()))

getfailed(x::ParallelJobs) = ParallelJobs(_selectby(x, Failed()))

function _selectby(j::ParallelJobs, st::JobStatus)  # Do not export!
    res = AtomicJob[]
    for subjob in j.subjobs
        if getstatus(subjob) === st
            push!(res, subjob)
        end
    end
    return res
end # function _selectby

function Base.show(io::IO, x::ParallelJobs)
    n = length(x.subjobs)
    println(io, "# $n subjobs in this job:")
    for (i, subjob) in enumerate(x.subjobs)
        print(io, lpad("[$i", ndigits(n) + 2), "] ", _emoji(subjob))
        printstyled(io, " ", subjob.cmd; bold = true)
        printstyled(
            io,
            " @ ",
            format(subjob.starttime, "Y/mm/dd H:M:S"),
            ", uses ",
            _readabletime(timecost(subjob)),
            '\n';
            color = :light_black,
        )
    end
end # function Base.show

_emoji(subjob) = isrunning(subjob) ? '🚧' : issucceeded(subjob) ? '✅' : '❌'

_readabletime(t) = canonicalize(CompoundPeriod(t))  # Do not export!

function div_nprocs(np, nj)
    quotient, remainder = divrem(np, nj)
    if !iszero(remainder)
        @warn "The processes are not fully balanced! Consider the number of subjobs!"
    end
    return quotient
end # function div_nprocs

function peepstdout(x::AtomicJob)
    out = getstdout(x.cmd)
    if out !== nothing
        println(read(out, String))
    end
end # function peepstdout
peepstdout(x::ParallelJobs) = foreach(peepstdout, x.subjobs)

function peepstdin(x::AtomicJob)
    out = getstdin(x.cmd)
    if out !== nothing
        println(read(out, String))
    end
end # function peepstdin
peepstdin(x::ParallelJobs) = foreach(peepstdin, x.subjobs)

function peepstderr(x::AtomicJob)
    out = getstderr(x.cmd)
    if out !== nothing
        println(read(out, String))
    end
end # function peepstdin
peepstderr(x::ParallelJobs) = foreach(peepstderr, x.subjobs)

getstdin(x::Base.CmdRedirect) = x.stream_no == 0 ? x.handle.filename : getstdin(x.cmd)
getstdin(::Base.AbstractCmd) = nothing

getstdout(x::Base.CmdRedirect) = x.stream_no == 1 ? x.handle.filename : getstdout(x.cmd)
getstdout(::Base.AbstractCmd) = nothing

getstderr(x::Base.CmdRedirect) = x.stream_no == 2 ? x.handle.filename : getstderr(x.cmd)
getstderr(::Base.AbstractCmd) = nothing

Base.iterate(x::ParallelJobs) = iterate(x.subjobs)
Base.iterate(x::ParallelJobs, state) = iterate(x.subjobs, state)

Base.getindex(x::ParallelJobs, i) = getindex(x.subjobs, i)

Base.firstindex(x::ParallelJobs) = 1

Base.lastindex(x::ParallelJobs) = length(x.subjobs)

end
