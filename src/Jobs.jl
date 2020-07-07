module Jobs

using Dates: DateTime, CompoundPeriod, now, canonicalize
using Distributed

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
    getfailed

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
function launchjob(tracker::JobTracker, interval = 3)
    return JobTracker(map(tracker.subjobs) do subjob
        if getstatus(subjob) ∈ (RUNNING, SUCCEEDED)
            subjob
        else
            sleep(interval)
            x = OneShot(subjob.cmd)
            x.status = RUNNING
            x.ref = @spawn begin
                x.starttime = now()
                try
                    z = run(subjob.cmd; wait = true)
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
    end)
end # function launchjob

getstatus(x::OneShot) = x.status
getstatus(x::JobTracker) = map(getstatus, x.subjobs)

starttime(x::OneShot) = x.starttime
starttime(x::JobTracker) = map(starttime, x.subjobs)

stoptime(x::OneShot) = isrunning(x) ? nothing : x.stoptime
stoptime(x::JobTracker) = map(stoptime, x.subjobs)

timecost(x::OneShot) = isrunning(x) ? nothing : x.stoptime - x.starttime
timecost(x::JobTracker) = map(timecost, x.subjobs)

getresult(x::OneShot) = isrunning(x) ? nothing : fetch(x.ref)
getresult(x::JobTracker) = map(getresult, x.subjobs)

isrunning(x::OneShot) = getstatus(x) === RUNNING

issucceeded(x::OneShot) = getstatus(x) === SUCCEEDED

isfailed(x::OneShot) = getstatus(x) === FAILED

getrunning(x::JobTracker) = JobTracker(_selectby(x, RUNNING))

getsucceeded(x::JobTracker) = JobTracker(_selectby(x, SUCCEEDED))

getfailed(x::JobTracker) = JobTracker(_selectby(x, FAILED))

function _selectby(j::JobTracker, st::JobStatus)  # Do not export!
    res = OneShot[]
    for subjob in j.subjobs
        if getstatus(subjob) === st
            push!(res, subjob)
        end
    end
    return res
end # function _selectby

function Base.show(io::IO, x::JobTracker)
    n = length(x.subjobs)
    println(io, "# $n subjobs in this job:")
    for (i, subjob) in enumerate(x.subjobs)
        println(io, lpad("[$i", ndigits(n) + 1), "] ", subjob.cmd)
        print(io, ' '^(ndigits(n) + 4))
        if isrunning(subjob)
            printstyled(
                io,
                "running for: ",
                _readabletime(now() - subjob.starttime),
                '\n';
                color = :blue,
            )  # Blue text
        elseif isfailed(subjob)
            printstyled(
                io,
                "failed after: ",
                _readabletime(subjob.stoptime - subjob.starttime),
                '\n';
                color = :red,
            )  # Red text
        elseif issucceeded(subjob)
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
