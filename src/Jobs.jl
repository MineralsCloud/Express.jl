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
    return JobTracker(pairs(tasks), [], [])
end # function launchjob
# function launchjob(cmds, environment::DockerEnvironment)
#     return map(cmds) do cmd
#         @async exec_run(environment.container, cmd; demux = true)
#     end
# end # function launchjob

function update!(x::JobTracker)
    @assert isvalid(x)
    for (i, task) in enumerate(x.running)
        if isready(task)
            try
                result = fetch(task)
                push!(x.succeeded, i => result)
            catch e
                push!(x.failed, i => e)
            finally
                deleteat!(x.running, i)
            end
        end
    end
    return x
end # function update!

Base.isvalid(x::JobTracker) =
    x.n == length(x.running) + length(x.succeeded) + length(x.failed)

end
