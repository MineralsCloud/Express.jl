module Jobs

using Dates: DateTime, CompoundPeriod, now, canonicalize, format
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
    getfailed,
    peepstdin,
    peepstdout,
    peepstderr

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

timecost(x::OneShot) = isrunning(x) ? now() - x.starttime : x.stoptime - x.starttime
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
        print(io, lpad("[$i", ndigits(n) + 2), "] ", _emoji(subjob))
        printstyled(io, " ", subjob.cmd; bold = true)
        printstyled(
            io,
            " @ ",
            format(subjob.starttime, "Y/mm/dd H:M:S"),
            ", uses ",
            timecost(subjob),
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

function peepstdout(x::OneShot)
    out = getstdout(x.cmd)
    if out !== nothing
        println(read(out, String))
    end
end # function peepstdout
peepstdout(x::JobTracker) = foreach(peepstdout, x.subjobs)

function peepstdin(x::OneShot)
    out = getstdin(x.cmd)
    if out !== nothing
        println(read(out, String))
    end
end # function peepstdin
peepstdin(x::JobTracker) = foreach(peepstdin, x.subjobs)

function peepstderr(x::OneShot)
    out = getstderr(x.cmd)
    if out !== nothing
        println(read(out, String))
    end
end # function peepstdin
peepstderr(x::JobTracker) = foreach(peepstderr, x.subjobs)

getstdin(x::Base.CmdRedirect) = x.stream_no == 0 ? x.handle.filename : getstdin(x.cmd)
getstdin(::Base.AbstractCmd) = nothing

getstdout(x::Base.CmdRedirect) = x.stream_no == 1 ? x.handle.filename : getstdout(x.cmd)
getstdout(::Base.AbstractCmd) = nothing

getstderr(x::Base.CmdRedirect) = x.stream_no == 2 ? x.handle.filename : getstderr(x.cmd)
getstderr(::Base.AbstractCmd) = nothing

end
