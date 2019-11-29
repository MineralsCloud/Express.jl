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

function distribute_process(manager::ClusterManager, job::Job, worker_ids::Integer = workers())
    # mpirun -np $n pw.x -in $in -out $out
    addprocs(manager, job.directives...)
    futures = Future[]
    for id in worker_ids
        push!(futures, @spawnat(id, run(job.action)))
    end
    return futures
end # function submit

end
