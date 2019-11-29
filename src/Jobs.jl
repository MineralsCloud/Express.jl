module Jobs

export Job, SubJob, distribute_process, isjobdone

using Dates: DateTime, now
using Distributed

using ClusterManagers
using Parameters: @with_kw

@with_kw struct Job
    id::String
    action::Cmd
    name::String = ""
    priority::Int = 0
    time_created::DateTime = now()
    directives::Dict = Dict{String,Any}()
end

struct SubJob
    parent::Job
    worker_id::Int
    ref::Future
end

function distribute_process(manager::ClusterManager, job::Job, worker_ids = workers())
    # mpirun -np $n pw.x -in $in -out $out
    # Similar to `invoke_on_workers` in https://cosx.org/2017/08/distributed-learning-in-julia
    subjobs = Vector{SubJob}(undef, length(worker_ids))
    for (i, id) in enumerate(worker_ids)
        subjobs[i] = SubJob(job, id, @spawnat id run(job.action))
    end
    return subjobs
end # function distribute_process

function isjobdone(subjobs::AbstractArray{SubJob})
    return all(map(isready, subjobs))
end # function isjobdone

end
