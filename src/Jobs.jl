module Jobs

export Job,
    distribute_process

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

function distribute_process(manager::ClusterManager, job::Job, worker_ids = workers())
    # mpirun -np $n pw.x -in $in -out $out
    # Similar to `invoke_on_workers` in https://cosx.org/2017/08/distributed-learning-in-julia
    futures = Vector{Future}(undef, length(worker_ids))
    for (i, id) in enumerate(worker_ids)
        @async begin
            futures[i] = @spawnat id run(job.action)
        end
    end
    return futures
end # function distribute_process

end
