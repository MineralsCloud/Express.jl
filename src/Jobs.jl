module Jobs

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

function submit(manager::ClusterManager, job::Job)
    addprocs(manager, job.directives...)
    job_ids = Int[]
    for i in workers()
        push!(job_ids, fetch(@spawnat i run(``)))
    end
    for i in workers()
        rmprocs(i)
    end
end # function submit

end
