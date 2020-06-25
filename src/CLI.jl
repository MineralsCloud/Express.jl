module CLI

export MpiCmd

struct MpiCmd
    n::UInt
    bin
    host
    arch
    wdir
    path
    file
    configfile
    env
end
MpiCmd(
    n;
    bin = "mpiexec",
    host = "",
    arch = "",
    wdir = "",
    path = String[],
    file = "",
    configfile = "",
    env = ENV,
) = MpiCmd(n, bin, host, arch, wdir, path, file, configfile, env)
# The docs are from https://www.mpich.org/static/docs/v3.3/www1/mpiexec.html.

end # module CLI
