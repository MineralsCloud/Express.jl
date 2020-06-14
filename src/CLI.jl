module CLI

export mpicmd

function mpicmd(
    n,
    cmd;
    bin = "mpiexec",
    host = "",
    arch = "",
    wdir = "",
    path = String[],
    file = "",
    configfile = "",
    env = ENV,
)
    options = String[]
    for (f, v) in zip((:host, :arch, :wdir, :file, :configfile), (host, arch, wdir, file, configfile))
        if !isempty(v)
            push!(options, "-$f", v)
        end
    end
    return Cmd([bin, "-n", string(n), options..., cmd.exec...])
end # function mpicmd
# The docs are from https://www.mpich.org/static/docs/v3.3/www1/mpiexec.html.

end # module CLI
