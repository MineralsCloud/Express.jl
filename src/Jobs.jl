module Jobs

using Distributed

using ClusterManagers
using Parameters: @with_kw
using QuantumESPRESSOBase.CLI: PWCmd
using Setfield: @set!

export MpiExec, JobStatus
export nprocs_task,
    distribute_process, isjobdone, tasks_running, tasks_exited, fetch_results, jobstatus

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

@enum JobStatus begin
    PENDING
    RUNNING
    EXITED
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
        if isready(task)
            try
                ref = fetch(task)
                status[i] = success(ref) ? JobStatus(:succeeded) : JobStatus(:failed)
            catch e
                status[i] = JobStatus(:failed)
            end
        end
        status[i] = JobStatus(:running)
    end
    return ids, status
end # function jobstatus

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
