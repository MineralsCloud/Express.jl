module Jobs

using Dates: DateTime, CompoundPeriod, now, canonicalize
using Distributed

export div_nprocs, launchjob, running, succeeded, failed, starttime, endtime, usetime

struct JobTracker
    subjobs::Vector{Future}
    starttime::Vector{DateTime}
    endtime::Vector{DateTime}
end

function launchjob(cmds, interval = 3)
    tracker = JobTracker(
        vec(similar(cmds, Future)),
        vec(similar(cmds, DateTime)),
        vec(similar(cmds, DateTime)),
    )
    for (i, cmd) in enumerate(cmds)
        sleep(interval)
        tracker.subjobs[i] = @spawn begin
            tracker.starttime[i] = now()
            try
                x = run(cmd; wait = true)
            finally
                tracker.endtime[i] = now()
                x
            end
        end
    end
    return tracker
end # function launchjob

running(x::JobTracker) = filter(!isready, x.subjobs)

succeeded(x::JobTracker) =
    filter(y -> !isa(y, RemoteException), map(fetch, filter(isready, x.subjobs)))

failed(x::JobTracker) =
    filter(y -> y isa RemoteException, map(fetch, filter(isready, x.subjobs)))

starttime(x::JobTracker) = x.starttime

function endtime(x::JobTracker)
    return map(x.subjobs, x.endtime) do subjob, endtime
        if isready(subjob)
            endtime
        else
            nothing
        end
    end
end # function endtime

function usetime(x::JobTracker)
    return map(x.subjobs, x.endtime, x.starttime) do subjob, endtime, starttime
        if isready(subjob)
            endtime - starttime
        else
            nothing
        end
    end
end # function usetime

function Base.show(io::IO, x::JobTracker)
    n = length(x.subjobs)
    println(io, "# $n subjobs in this job:")
    for (i, (starttime, endtime, subjob)) in
        enumerate(zip(x.starttime, x.endtime, x.subjobs))
        println(io, " [", lpad(i, ndigits(n)), "] ", subjob)
        print(io, ' '^(ndigits(n) + 4))
        if !isready(subjob)
            printstyled(
                io,
                "running for: ",
                _readabletime(now() - starttime),
                '\n';
                color = :blue,
            )  # Blue text
        else
            res = fetch(subjob)
            if res isa RemoteException
                printstyled(
                    io,
                    "failed after: ",
                    _readabletime(endtime - starttime),
                    '\n';
                    color = :red,
                )  # Red text
            else
                printstyled(
                    io,
                    "succeeded after: ",
                    _readabletime(endtime - starttime),
                    '\n';
                    color = :green,
                )  # Green text
            end
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

end
