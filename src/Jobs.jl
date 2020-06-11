module Jobs

using Distributed
using DockerPy.Containers

export JobStatus, JobResult
export nprocs_task,
    distribute_process,
    isjobdone,
    tasks_running,
    tasks_exited,
    fetch_results,
    jobstatus,
    jobresult

@enum JobStatus begin
    RUNNING
    EXITED
end

@enum JobResult begin
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

function distribute_process(
    cmds::AbstractArray,
    ids::AbstractArray{<:Integer} = workers();
    isdocker::Bool = true,
    container,
    inputs,
)
    # Similar to `invoke_on_workers` in https://cosx.org/2017/08/distributed-learning-in-julia
    return map(cmds, ids, inputs) do cmd, id, input  # promises
        if isdocker
            exec_run(container, cmd; demux = true, workdir = dirname(input))
        else
            @spawnat id run(cmd, wait = true)  # TODO: Must wait?
        end
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

function jobstatus(bag::AbstractVector{Future})
    return map(bag) do task
        task.where => isready(task) ? EXITED : RUNNING  # id => status
    end
end # function jobstatus

function jobresult(bag::AbstractVector{Future})
    return map(bag) do task
        id = task.where
        if isready(task)
            try
                ref = fetch(task)
                result = success(ref) ? SUCCEEDED : FAILED
            catch e
                result = FAILED
            end
        else
            result = nothing
        end
        id => result
    end
end # function jobresult

function fetch_results(bag::AbstractVector)
    return map(bag) do task
        id = task.where
        if isready(task)
            try
                result = fetch(task)
            catch e
                result = e
            end
        else
            result = nothing
        end
        id => result
    end
end # function fetch_results

end
