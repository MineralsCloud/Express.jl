module Jobs

using Dates: DateTime, CompoundPeriod, now, canonicalize
using Distributed
using DockerPy.Containers

export div_nprocs, launchjob, running, succeeded, failed, starttime, endtime, usetime

struct JobTracker
    subjobs::Vector{Future}
    starttime::Vector{DateTime}
    endtime::Vector{DateTime}
end

function launchjob(cmds, interval = 3)
    tracker =
        JobTracker(similar(cmds, Future), similar(cmds, DateTime), similar(cmds, DateTime))
    for (i, cmd) in enumerate(cmds)
        sleep(interval)
        tracker.subjobs[i] = @spawn begin
            tracker.starttime[i] = now()
            x = run(cmd; wait = true)
            tracker.endtime[i] = now()
            x
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
    print(
        io,
        join(
            (
                "\033[34mrunning:\033[0m ",  # Blue text
                running(x),
                "\033[32msucceeded:\033[0m ",  # Green text
                succeeded(x),
                "\033[31mfailed:\033[0m ",  # Red text
                failed(x),
                "\033[34mtime used:\033[0m ",  # Green text
                map(
                    y -> y !== nothing ? canonicalize(CompoundPeriod(y)) : nothing,
                    usetime(x),
                ),
            ),
            '\n',
        ),
    )
end # function Base.show

function div_nprocs(np, nj)
    quotient, remainder = divrem(np, nj)
    if !iszero(remainder)
        @warn "The processes are not fully balanced! Consider the number of subjobs!"
    end
    return quotient
end # function div_nprocs

end
