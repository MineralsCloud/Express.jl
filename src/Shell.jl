module Shell

using SimpleWorkflow: InternalAtomicJob, Script, parallel

export @intjob

function distprocs(nprocs, njobs)
    quotient, remainder = divrem(nprocs, njobs)
    if !iszero(remainder)
        @warn "The processes are not fully balanced! Consider the number of subjobs!"
    end
    return quotient
end

macro intjob(ex)  # See https://github.com/JuliaLang/julia/blob/ab5853f/base/task.jl#L111-L113
    return :(InternalAtomicJob(() -> $(esc(ex))))
end

end
