module Jobs

using Distributed
using Setfield: @set!

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

function nprocs_task(total_num::Int, nsubjob::Int)
    quotient, remainder = divrem(total_num, nsubjob)
    if remainder != 0
        @warn("The processes are not fully balanced! Consider the number of subjobs!")
    end
    return quotient
end # function nprocs_task

function distribute_process(cmds::AbstractArray, ids::AbstractArray{<:Integer} = workers())
    # Similar to `invoke_on_workers` in https://cosx.org/2017/08/distributed-learning-in-julia
    if length(cmds) != length(ids)  # The size of them can be different, but not length.
        throw(DimensionMismatch("`cmds` has different length than `ids`!"))
    end
    promises = similar(cmds, Future)  # It can be of different size than `ids`!
    for (i, (cmd, id)) in enumerate(zip(cmds, ids))
        promises[i] = @spawnat id run(cmd, wait = true)  # TODO: Must wait?
    end
    return promises
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
    ids, status = Vector{Int}(undef, length(bag)), Vector{JobStatus}(undef, length(bag))
    for (i, task) in enumerate(bag)
        ids[i] = task.where
        status[i] = isready(task) ? EXITED : RUNNING
        if isready(task)
            status[i] = EXITED
        end
        status[i] = RUNNING
    end
    return ids, status
end # function jobstatus

function jobresult(bag::AbstractVector{Future})
    ids, results = Vector{Int}(undef, length(bag)),
    Vector{Union{JobResult,Nothing}}(undef, length(bag))
    for (i, task) in enumerate(bag)
        ids[i] = task.where
        results[i] = isready(task) ? EXITED : RUNNING
        if isready(task)
            try
                ref = fetch(task)
                results[i] = success(ref) ? SUCCEEDED : FAILED
            catch e
                results[i] = FAILED
            end
        end
        results[i] = nothing
    end
    return ids, results
end # function jobresult

function fetch_results(bag::AbstractVector)
    ids, results = Vector{Int}(undef, length(bag)), Vector{Any}(undef, length(bag))
    for (i, task) in enumerate(bag)
        ids[i] = task.where
        if isready(task)
            try
                results[i] = fetch(task)
            catch e
                results[i] = e
            end
        end
        results[i] = nothing
    end
    return ids, results
end # function fetch_results

end
