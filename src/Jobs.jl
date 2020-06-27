module Jobs

using Distributed
using DockerPy.Containers

export nprocs_task, launchjob, update!

struct JobTracker
    running
    succeeded
    failed
    n::Int
    JobTracker(running, succeeded, failed) = new(
        running,
        succeeded,
        failed,
        length(running) + length(succeeded) + length(failed),
    )
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
    @assert isvalid(x)
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
    foreach(
        x -> println(io, x),
        ("running:", x.running, "succeeded:", x.succeeded, "failed:", x.failed),
    )
end # function Base.show

Base.isvalid(x::JobTracker) =
    x.n == length(x.running) + length(x.succeeded) + length(x.failed)

end
