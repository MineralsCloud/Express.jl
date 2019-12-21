module Jobs

using Distributed

using ClusterManagers
using Parameters: @with_kw
using QuantumESPRESSOBase.CLI: PWCmd
using Setfield: @set!

export MpiExec, BagOfTasks, TaskStatus
export nprocs_task, distribute_process, isjobdone, fetch_results, jobstatus

@with_kw struct MpiExec <: Base.AbstractCmd
    # The docs are from https://www.mpich.org/static/docs/v3.3/www1/mpiexec.html.
    "The path to the executable, defaults to \"mpiexec\""
    which::String = "mpiexec"
    "Specify the number of processes to use"
    n::Int
    "Name of host on which to run processes"
    host::String = ""
    "Pick hosts with this architecture type"
    arch::String = ""
    "`cd` to this one before running executable"
    wdir::String = ""
    "Use this to find the executable"
    path::Vector{String} = []
    "Implementation-defined specification file"
    file::String = ""
    configfile::String = ""
    cmd::Base.AbstractCmd
    "Set environment variables to use when running the command, defaults to `ENV`"
    env = ENV
end

struct TaskStatus{T}
    function TaskStatus{T}() where {T}
        T ∈ (:pending, :running, :succeeded, :failed) ||
        throw(ArgumentError("the task status should be one of `:pending`, `:running`, `:succeeded`, `:failed`!"))
        return new()
    end # function TaskStatus
end
TaskStatus(T) = TaskStatus{T}()

struct BagOfTasks{T<:AbstractArray}
    tasks::T
end

function nprocs_task(total_num::Int, nsubjob::Int)
    quotient, remainder = divrem(total_num, nsubjob)
    if remainder != 0
        @warn("The processes are not fully balanced! Consider the number of subjobs!")
    end
    return quotient
end # function nprocs_task

function distribute_process(
    cmds::AbstractArray{T},
    ids::AbstractArray{<:Integer} = workers(),
) where {T<:MpiExec}
    # Similar to `invoke_on_workers` in https://cosx.org/2017/08/distributed-learning-in-julia
    if length(cmds) != length(ids)  # The size of them can be different, but not length.
        throw(DimensionMismatch("`cmds` has different length than `ids`!"))
    end
    promises = similar(cmds, Future)  # It can be of different size than `ids`!
    for (i, (cmd, id)) in enumerate(zip(cmds, ids))
        promises[i] = @spawnat id run(convert(Cmd, cmd), wait = true)  # TODO: Must wait?
    end
    return BagOfTasks(promises)
end # function distribute_process

function isjobdone(bag::BagOfTasks)
    return all(map(isready, bag.tasks))
end # function isjobdone

function tasks_running(bag::BagOfTasks)
    return filter(!isready, bag.tasks)
end # function tasks_running

function tasks_exited(bag::BagOfTasks)
    return map(fetch, filter(isready, bag.tasks))
end # function subjobs_exited

function jobstatus(bag::BagOfTasks)
    map(bag.tasks) do task
        if isready(task)
            return success(fetch(task)) ? TaskStatus(:succeeded) : TaskStatus(:failed)
        end
        return isempty(task) ? TaskStatus(:pending) : TaskStatus(:running)
    end
end # function jobstatus

function fetch_results(bag::BagOfTasks)
    return map(bag) do x
        if isready(x)
            try
                fetch(x)
            catch e
                e
            end
        end
    end
end # function fetch_results

function Base.convert(::Type{Cmd}, cmd::MpiExec)
    options = String[]
    # for f in fieldnames(typeof(cmd))[3:end]  # Join options
    #     v = getfield(cmd, f)
    #     if !iszero(v)
    #         push!(options, string(" -", f, ' ', v))
    #     else
    #         push!(options, "")
    #     end
    # end
    return Cmd(
        `$(cmd.which) -np $(cmd.n) $(options...) $(convert(Cmd, cmd.cmd))`,
        env = cmd.env,
        dir = cmd.wdir,
    )
end # function Base.convert

end
