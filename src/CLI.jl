module CLI

using Parameters: @with_kw_noshow
using QuantumESPRESSO.CLI: PWCmd
using REPL.Terminals: AbstractTerminal

export MpiExec, DockerExec

@with_kw_noshow struct DockerExec <: Base.AbstractCmd
    which::String = "docker"
    container::String
    cmd::Base.AbstractCmd
    detach::Bool = false
    env::Base.EnvDict = ENV
    interactive::Bool = false
    tty::IO = stdin
    user::UInt64 = 0
    workdir::String = pwd()
end

@with_kw_noshow struct MpiExec <: Base.AbstractCmd
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
    wdir::String = pwd()
    "Use this to find the executable"
    path::Vector{String} = []
    "Implementation-defined specification file"
    file::String = ""
    configfile::String = ""
    cmd::Base.AbstractCmd
    "Set environment variables to use when running the command, defaults to `ENV`"
    env = ENV
end

function Base.convert(::Type{Cmd}, cmd::MpiExec)
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
        `$(cmd.which) -np $(cmd.n) $(options...) $(convert(Cmd, cmd.cmd).exec)`,
        env = cmd.env,
        dir = cmd.wdir,
    )
end # function Base.convert
function Base.convert(::Type{Cmd}, cmd::DockerExec)
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
        `$(cmd.which) exec $(options...) $(cmd.container) $(convert(Cmd, cmd.cmd).exec)`,
        env = cmd.env,
        dir = cmd.workdir,
    )
end # function Base.convert

end # module CLI

