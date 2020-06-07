module Jobs

using Distributed
using DockerPy.Containers

using ..Environments: DockerEnvironment, LocalEnvironment

export JobStatus
export nprocs_task,
    distribute_process,
    isjobdone,
    tasks_running,
    tasks_exited,
    jobresult,
    jobstatus

@enum JobStatus begin
    RUNNING
    SUCCEEDED
    FAILED
end

function nprocs_task(total_num, nsubjob)
    quotient, remainder = divrem(total_num, nsubjob)
    if !iszero(remainder)
        @warn "The processes are not fully balanced! Consider the number of subjobs!"
    end
    return quotient
end # function nprocs_task

function distribute_process(cmds, ::LocalEnvironment)
    # Similar to `invoke_on_workers` in https://cosx.org/2017/08/distributed-learning-in-julia
    ids = addprocs(length(cmds))
    return map(cmds, ids) do cmd, id  # promises
        @spawnat id run(cmd; wait = true)  # Must wait, or else lose I/O streams
    end
end # function distribute_process
function distribute_process(cmds, environment::DockerEnvironment)
    # Similar to `invoke_on_workers` in https://cosx.org/2017/08/distributed-learning-in-julia
    ids = addprocs(length(cmds))
    return map(cmds, ids) do cmd, id  # promises
        exec_run(environment.container, cmd; demux = true)
    end
end # function distribute_process

function isjobdone(bag::AbstractVector)
    return all(map(isready, bag))
end # function isjobdone

function tasks_running(bag::AbstractVector)
    return filter(!isready, bag)
end # function tasks_running

function tasks_exited(bag::AbstractVector)
    return filter(isready, bag)
end # function subjobs_exited

function jobstatus(job::Future)
    id = job.where
    if isready(job)
        try
            ref = fetch(job)
            status = success(ref) ? SUCCEEDED : FAILED
        catch e
            status = FAILED
        end
    else
        status = RUNNING
    end
    return id => status
end # function jobstatus
jobstatus(bag) = map(jobstatus, bag)

function jobresult(job::Future)
    id = job.where
    if isready(job)
        try
            result = fetch(job)
        catch e
            result = e
        end
    else
        result = nothing
    end
    return id => result
end # function jobresult
jobresult(bag) = map(jobresult, bag)

end
