module Jobs

using Compat: isnothing
using Distributed
using DockerPy.Containers

using ..Environments: DockerEnvironment, LocalEnvironment

export JobStatus
export nprocs_task, distributejobs, isjobdone, jobresult, jobstatus

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

function distributejobs(cmds, ::LocalEnvironment)
    return map(cmds) do cmd
        @async run(cmd; wait = true)  # Must wait, or else lose I/O streams
    end
end # function distributejobs
function distributejobs(cmds, environment::DockerEnvironment)
    return map(cmds) do cmd
        @async exec_run(environment.container, cmd; demux = true)
    end
end # function distributejobs

isjobdone(job::Future) = isready(job)
isjobdone(bag) = map(isjobdone, bag)

function jobstatus(job::Future)
    id, result = jobresult(job)
    return isnothing(result) ? RUNNING : success(result) ? SUCCEEDED : FAILED
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
