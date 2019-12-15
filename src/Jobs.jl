module Jobs

using Distributed

using ClusterManagers
using Parameters: @with_kw
using Setfield: @set!

using QuantumESPRESSOBase.CLI

using Express

export MpiExec
export nprocs_per_subjob, distribute_process, isjobdone, fetch_results

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
    subcmd::Base.AbstractCmd
    "Set environment variables to use when running the command, defaults to `ENV`"
    env::Base.EnvDict = ENV  # FIXME: What is this type?
end

struct BagOfTasks{T<:AbstractArray}
    tasks::T
end

function nprocs_per_subjob(total_num::Int, nsubjob::Int)
    quotient, remainder = divrem(total_num, nsubjob)
    if remainder != 0
        @warn("The processes are not fully balanced! Consider the number of subjobs!")
    end
    return quotient
end # function nprocs_per_subjob

function distribute_process(
    cmds::AbstractArray{T},
    ids::AbstractArray{<:Integer} = workers(),
) where {T<:MpiExec}
    # Similar to `invoke_on_workers` in https://cosx.org/2017/08/distributed-learning-in-julia
    if size(cmds) != size(ids)
        throw(DimensionMismatch("`cmds` has different size than `ids`!"))
    end
    refs = similar(ids, Future)
    for (i, (cmd, id)) in enumerate(zip(cmds, ids))
        refs[i] = @spawnat id run(Cmd(cmd), wait = true)  # TODO: Must wait?
    end
    return refs
end # function distribute_process

function isjobdone(refs::AbstractArray{Future})
    return all(map(isready, refs))
end # function isjobdone

function subjobs_running(refs::AbstractArray{Future})
    return filter(!isready, refs)
end # function monitor

function subjobs_exited(refs::AbstractArray{Future})
    return map(fetch, filter(isready, refs))
end # function subjobs_exited

function fetch_results(refs::AbstractArray{Future})
    return map(refs) do x
        if isready(x)
            try
                fetch(x)
            catch e
                e
            end
        end
    end
end # function fetch_results

Base.Cmd(cmd::MpiExec) = Cmd(`$(cmd.which) -np $(cmd.n) $(Cmd(cmd.subcmd))`, env = ENV, dir = cmd.wdir)

end
