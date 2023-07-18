export procs_per_job

function procs_per_job(nprocs, jobsize)
    quotient, remainder = divrem(nprocs, jobsize)
    if !iszero(remainder)
        @warn "the processes are not fully balanced! Consider the number of subjobs!"
    end
    return quotient
end
