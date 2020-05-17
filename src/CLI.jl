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
    wdir = "",
    path = [],
    file = "",
    configfile = "",
    env = ENV,
) = MpiExec(n, which, host, arch, wdir, path, file, configfile, env)
function (exec::MpiExec)(cmd::Base.AbstractCmd)
    options = String[]
    for f in (:host, :arch, :wdir, :file, :configfile)
        v = getfield(exec, f)
        if !isempty(v)
            push!(options, "-$f", v)
        end
    end
    _deepinject!(cmd, [exec.which, "-np", string(exec.n), options...])  # Horrible
    # return Cmd(cmd, env = exec.env, dir = exec.wdir)
    return cmd
end

_deepinject!(cmd::Cmd, prefix::AbstractVector) = prepend!(cmd.exec, prefix)
_deepinject!(cmd::Base.CmdRedirect, prefix::AbstractVector) = _deepinject!(cmd.cmd, prefix)

end # module CLI
