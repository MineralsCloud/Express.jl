module Jobs

using Dates: DateTime, now
using Distributed

using ClusterManagers
using Parameters: @with_kw
using Setfield: @set!

using QuantumESPRESSOBase.CLI

using Express

export MpiCmd, Job
export nprocs_per_subjob, distribute_process, isjobdone, fetch_results

@with_kw struct MpiCmd
    exec::String = "mpirun"
    np::Int
    subcmd::CLI.QuantumESPRESSOCmd
end

@with_kw struct Job
    id::String
    action::Cmd
    name::String = ""
    priority::Int = 0
    time_created::DateTime = now()
    directives::Dict = Dict{String,Any}()
end

function nprocs_per_subjob(total_num::Int, nsubjob::Int)
    quotient, remainder = divrem(total_num, nsubjob)
    if remainder != 0
        @warn("The processes are not fully balanced! Consider the number of subjobs!")
    end
    return quotient
end # function nprocs_per_subjob

function distribute_process(cmd::MpiCmd, worker_ids = workers())
    # mpirun -np $n pw.x -in $in -out $out
    # Similar to `invoke_on_workers` in https://cosx.org/2017/08/distributed-learning-in-julia
    refs = Vector{Future}(undef, nworkers())
    @set! cmd.np = nprocs_per_subjob(cmd.np, length(worker_ids))
    for (i, id) in enumerate(worker_ids)
        refs[i] = @spawnat id run(commandify(cmd))
    end
    return refs
end # function distribute_process

function isjobdone(refs::AbstractArray{Future})
    return all(map(x -> isready(x), refs))
end # function isjobdone

function fetch_results(refs::AbstractArray{Future})
    return map(refs) do x
        isready(x) ? fetch(x) : nothing
    end
end # function fetch_results

CLI.commandify(cmd::MpiCmd) = `$(cmd.exec) -np $(cmd.np) $(commandify(cmd.subcmd))`

end
