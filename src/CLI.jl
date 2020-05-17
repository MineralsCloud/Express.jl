module CLI

export MpiExec, DockerExec

mutable struct DockerExec
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
    container;
    which = "docker",
    detach = false,
    env = ENV,
    interactive = false,
    tty = stdin,
    user = 0,
    workdir = pwd(),
) = DockerExec(container, which, detach, env, interactive, tty, user, workdir)
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
        Cmd([
            exec.which,
            "exec",
            options...,
            exec.container,
            "sh",
            "-c",
        ]),
        env = exec.env,
        dir = exec.workdir,
    )
end # function Base.Cmd
(exec::DockerExec)(cmd::Base.AbstractCmd) = Cmd([string(Cmd(exec)), "sh", "-c", string(cmd)[2:end-1]])

mutable struct MpiExec
    # The docs are from https://www.mpich.org/static/docs/v3.3/www1/mpiexec.html.
    "Specify the number of processes to use"
    n::Int
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
    env::Any
end
MpiExec(
    n;
    which = "mpiexec",
    host = "",
    arch = "",
    wdir = pwd(),
    path = [],
    file = "",
    configfile = "",
    env = ENV,
) = MpiExec(n, which, host, arch, wdir, path, file, configfile, env)
(exec::MpiExec)(cmd::Base.AbstractCmd) = `$(Cmd(exec)) $cmd`
(exec::MpiExec)(cmd::Base.CmdRedirect) = pipeline(`$(Cmd(exec)) $(cmd.cmd)`)

function Base.Cmd(exec::MpiExec)
    options = String[]
    # for f in Iterators.drop(fieldnames(typeof(exec)), 3)  # 4 to end
    #     v = getfield(exec, f)
    #     if !iszero(v)
    #         push!(options, string(" -", f, ' ', v))
    #     else
    #         push!(options, "")
    #     end
    # end
    return Cmd(
        Cmd([
            exec.which,
            "-np",
            string(exec.n),
            options...,
        ]),
        env = exec.env,
        dir = exec.wdir,
    )
end # function Base.Cmd


end # module CLI
