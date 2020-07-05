module Jobs

using Dates: DateTime, CompoundPeriod, now, canonicalize
using Distributed
using DockerPy.Containers

export nprocs_task, launchjob, update!

struct JobTracker
    running
    succeeded
    failed
    starttime::DateTime
    JobTracker(running, succeeded, failed) = new(running, succeeded, failed, now())
end

function nprocs_task(total_num, nsubjob)
    quotient, remainder = divrem(total_num, nsubjob)
    if !iszero(remainder)
        @warn "The processes are not fully balanced! Consider the number of subjobs!"
    end
    return quotient
end # function nprocs_task

function launchjob(cmds; sleepfor = 5)
    tasks = map(cmds) do cmd
        sleep(sleepfor)
        @spawn run(cmd; wait = true)  # Must wait, or else lose I/O streams
    end
    return JobTracker(Dict(pairs(tasks)), Dict(), Dict())
end # function launchjob

function update!(x::JobTracker)
    for (i, task) in x.running
        if isready(task)
            result = fetch(task)
            if result isa RemoteException
                push!(x.failed, i => result)
            else
                push!(x.succeeded, i => result)
            end
            pop!(x.running, i)
        end
    end
    return x
end # function update!

function Base.show(io::IO, x::JobTracker)
    println(io, "time has passed: ", canonicalize(CompoundPeriod(now() - x.starttime)))
    foreach(
        x -> println(io, x),
        ("running:", x.running, "succeeded:", x.succeeded, "failed:", x.failed),
    )
end # function Base.show

end
