module Jobs

using Dates: DateTime, CompoundPeriod, now, canonicalize
using Distributed

export div_nprocs,
    launchjob,
    starttime,
    stoptime,
    usetime,
    status,
    running,
    succeeded,
    failed,
    getstdin,
    getstdout,
    getstderr

@enum JobStatus begin
    RUNNING
    SUCCEEDED
    FAILED
end

mutable struct OneShot
    cmd::Base.AbstractCmd
    ref::Future
    starttime::DateTime
    stoptime::DateTime
    status::JobStatus
    OneShot(cmd) = new(cmd)
end

struct JobTracker
    subjobs::Vector{OneShot}
end

function launchjob(cmds, interval = 3)
    subjobs = map(enumerate(cmds)) do (i, cmd)
        sleep(interval)
        x = OneShot(cmd)
        x.status = RUNNING
        x.ref = @spawn begin
            x.starttime = now()
            try
                z = run(cmd; wait = true)
                x.status = SUCCEEDED
            catch
                x.status = FAILED
            finally
                x.stoptime = now()
                z
            end
        end
        x
    end
    return JobTracker(vec(subjobs))
end # function launchjob

status(x::OneShot) = x.status
status(x::JobTracker) = map(status, x.subjobs)

starttime(x::OneShot) = x.starttime
starttime(x::JobTracker) = map(starttime, x.subjobs)

stoptime(x::OneShot) = status(x) == RUNNING ? nothing : x.stoptime
stoptime(x::JobTracker) = map(stoptime, x.subjobs)

usetime(x::OneShot) = status(x) == RUNNING ? nothing : x.stoptime - x.starttime
usetime(x::JobTracker) = map(usetime, x.subjobs)

running(x::JobTracker) =
    map(Base.Fix1(getindex, x.subjobs), findall(==(RUNNING), status(x)))

succeeded(x::JobTracker) =
    map(Base.Fix1(getindex, x.subjobs), findall(==(SUCCEEDED), status(x)))

failed(x::JobTracker) = map(Base.Fix1(getindex, x.subjobs), findall(==(FAILED), status(x)))

function Base.show(io::IO, x::JobTracker)
    n = length(x.subjobs)
    println(io, "# $n subjobs in this job:")
    for (i, subjob) in enumerate(x.subjobs)
        println(io, lpad("[$i", ndigits(n) + 1), "] ", subjob.cmd)
        print(io, ' '^(ndigits(n) + 4))
        status = status(subjob)
        if status == RUNNING
            printstyled(
                io,
                "running for: ",
                _readabletime(now() - subjob.starttime),
                '\n';
                color = :blue,
            )  # Blue text
        elseif status == FAILED
            printstyled(
                io,
                "failed after: ",
                _readabletime(subjob.stoptime - subjob.starttime),
                '\n';
                color = :red,
            )  # Red text
        elseif status == SUCCEEDED
            printstyled(
                io,
                "succeeded after: ",
                _readabletime(subjob.stoptime - subjob.starttime),
                '\n';
                color = :green,
            )  # Green text
        else
            error("unknown status!")  # This should never happen!
        end
    end
end # function Base.show

_readabletime(t) = canonicalize(CompoundPeriod(t))  # Do not export!

function div_nprocs(np, nj)
    quotient, remainder = divrem(np, nj)
    if !iszero(remainder)
        @warn "The processes are not fully balanced! Consider the number of subjobs!"
    end
    return quotient
end # function div_nprocs

getstdin(x::Base.CmdRedirect) = x.stream_no == 0 ? x.handle.filename : getstdin(x.cmd)
getstdin(::Base.AbstractCmd) = nothing

getstdout(x::Base.CmdRedirect) = x.stream_no == 1 ? x.handle.filename : getstdout(x.cmd)
getstdout(::Base.AbstractCmd) = nothing

getstderr(x::Base.CmdRedirect) = x.stream_no == 2 ? x.handle.filename : getstderr(x.cmd)
getstderr(::Base.AbstractCmd) = nothing

end
