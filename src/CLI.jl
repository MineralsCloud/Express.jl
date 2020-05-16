module CLI

using Parameters: @with_kw_noshow
using QuantumESPRESSO.CLI: PWExec
using REPL.Terminals: AbstractTerminal

export MpiExec, DockerExec

mutable struct DockerExec
    cmd
    container::String
    which::String
    detach::Bool
    env::Base.EnvDict
    interactive::Bool
    tty::IO
    user::UInt
    workdir::String
end
DockerExec(
    cmd,
    container;
    which = "docker",
    detach = false,
    env = ENV,
    interactive = false,
    tty = stdin,
    user = 0,
    workdir = pwd(),
) = DockerExec(cmd, container, which, detach, env, interactive, tty, user, workdir)

mutable struct MpiExec
    # The docs are from https://www.mpich.org/static/docs/v3.3/www1/mpiexec.html.
    "Specify the number of processes to use"
    n::Int
    "Specify the command to run"
    cmd
    "The path to the executable, defaults to \"mpiexec\""
    which::String
    "Name of host on which to run processes"
    host::String
    "Pick hosts with this architecture type"
    arch::String
    "`cd` to this one before running executable"
    wdir::String
    "Use this to find the executable"
    path::Vector{String}
    "Implementation-defined specification file"
    file::String
    configfile::String
    "Set environment variables to use when running the command, defaults to `ENV`"
    env
end
MpiExec(
    n,
    cmd;
    which = "mpiexec",
    host = "",
    arch = "",
    wdir = pwd(),
    path = [],
    file = "",
    configfile = "",
    env = ENV,
) = MpiExec(n, cmd, which, host, arch, wdir, path, file, configfile, env)

function Base.Cmd(exec::MpiExec)
    options = String[]
    # for f in fieldnames(typeof(cmd))[3:end]  # Join options
    #     v = getfield(cmd, f)
    #     if !iszero(v)
    #         push!(options, string(" -", f, ' ', v))
    #     else
    #         push!(options, "")
    #     end
    # end
    return Cmd(
        `$(exec.which) -np $(exec.n) $(options...) $(Cmd(exec.cmd))`,
        env = exec.env,
        dir = exec.wdir,
    )
end # function Base.Cmd
function Base.Cmd(exec::DockerExec)
    options = String[]
    # for f in fieldnames(typeof(cmd))[3:end]  # Join options
    #     v = getfield(cmd, f)
    #     if !iszero(v)
    #         push!(options, string(" -", f, ' ', v))
    #     else
    #         push!(options, "")
    #     end
    # end
    return Cmd(
        `$(exec.which) exec $(options...) $(exec.container) $(Cmd(exec.cmd))`,
        env = exec.env,
        dir = exec.workdir,
    )
end # function Base.Cmd

end # module CLI
