using ExpressBase: Action

export distribute_procs

function distribute_procs(nprocs, jobsize)
    quotient, remainder = divrem(nprocs, jobsize)
    if !iszero(remainder)
        @warn "the processes are not fully balanced! Consider the number of subjobs!"
    end
    return quotient
end

function jobify end

struct ExpandConfig{T} <: Action{T} end  # To be extended
