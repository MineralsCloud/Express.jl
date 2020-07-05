module Jobs

using Dates: DateTime, CompoundPeriod, now, canonicalize
using Distributed
using DockerPy.Containers

export div_nprocs, launchjob, running, succeeded, failed

struct JobTracker
    subjobs
    starttime::DateTime
    JobTracker(subjobs) = new(subjobs, now())
end

function div_nprocs(total_num, nsubjob)
    quotient, remainder = divrem(total_num, nsubjob)
    if !iszero(remainder)
        @warn "The processes are not fully balanced! Consider the number of subjobs!"
    end
    return quotient
end # function div_nprocs

function launchjob(cmds, interval = 3)
    return JobTracker(map(cmds) do cmd
        sleep(interval)
        @spawn run(cmd; wait = true)  # Must wait, or else lose I/O streams
    end)
end # function launchjob

running(x::JobTracker) = filter(!isready, x.subjobs)
succeeded(x::JobTracker) = filter(success, map(fetch, filter(isready, x.subjobs)))
failed(x::JobTracker) = filter(!success, map(fetch, filter(isready, x.subjobs)))

function Base.show(io::IO, x::JobTracker)
    print(
        io,
        join(
            (
                "\033[34mtime elapsed:\033[0m ",  # Green text
                canonicalize(CompoundPeriod(now() - x.starttime)),
                "\033[34mrunning:\033[0m ",
                running(x),
                "\033[32msucceeded:\033[0m ",
                succeeded(x),
                "\033[31mfailed:\033[0m ",
                failed(x),
            ),
            '\n',
        ),
    )
end # function Base.show

end
