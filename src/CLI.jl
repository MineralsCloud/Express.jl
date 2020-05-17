module CLI

export MpiExec

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
