module Jobs

using Distributed
using DockerPy.Containers

using ..Workspaces: DockerWorkspace, LocalWorkspace

export JobStatus
export nprocs_task,
    distribute_process,
    isjobdone,
    tasks_running,
    tasks_exited,
    fetch_results,
    jobstatus

@enum JobStatus begin
    PENDING
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

function distribute_process(cmds, ::LocalWorkspace)
    # Similar to `invoke_on_workers` in https://cosx.org/2017/08/distributed-learning-in-julia
    ids = addprocs(length(cmds))
    return map(cmds, ids) do cmd, id  # promises
        @spawnat id run(cmd; wait = true)  # TODO: Must wait?
    end
end # function distribute_process
function distribute_process(cmds, workspace::DockerWorkspace)
    # Similar to `invoke_on_workers` in https://cosx.org/2017/08/distributed-learning-in-julia
    ids = addprocs(length(cmds))
    return map(cmds, ids) do cmd, id  # promises
        exec_run(workspace.container, cmd; demux = true)
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
        id = task.where
        if isready(task)
            try
                ref = fetch(task)
                status = success(ref) ? SUCCEEDED : FAILED
            catch e
                status = FAILED
            end
        else
            status = RUNNING
        end
        id => status
    end
end # function jobstatus

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
