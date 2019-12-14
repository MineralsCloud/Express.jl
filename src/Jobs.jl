module Jobs

using Dates: DateTime, now
using Distributed

using ClusterManagers
using Parameters: @with_kw
using Setfield: @set!

using QuantumESPRESSOBase.CLI

using Express

export MpiCmd
export nprocs_per_subjob, distribute_process, isjobdone, fetch_results

@with_kw struct MpiCmd <: Base.AbstractCmd
    exec::String = "mpirun"
    np::Int
    subcmd::CLI.QuantumESPRESSOCmd
    env::Base.EnvDict = ENV
    stdin = nothing
    stdout = nothing
    stderr = nothing
    append::Bool = false
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
) where {T<:MpiCmd}
    # mpirun -np $n pw.x -in $in -out $out
    # Similar to `invoke_on_workers` in https://cosx.org/2017/08/distributed-learning-in-julia
    if size(cmds) != size(ids)
        throw(DimensionMismatch("`cmds` has different size than `ids`!"))
    end
    refs = similar(ids, Future)
    for (i, id) in enumerate(ids)
        refs[i] = @spawnat id run(Cmd(cmds[i]), wait = false)
    end
    return refs
end # function distribute_process

function isjobdone(refs::AbstractArray{Future})
    return all(map(isready, refs))
end # function isjobdone

function fetch_results(refs::AbstractArray{Future})
    return map(refs) do x
        isready(x) ? fetch(x) : nothing
    end
end # function fetch_results

Base.Cmd(cmd::MpiCmd) = pipeline(
    Cmd(`$(cmd.exec) -np $(cmd.np) $(Cmd(cmd.subcmd))`, env = ENV),
    cmd.stdin,
    cmd.stdout,
    cmd.stderr,
    cmd.append,
)

end
