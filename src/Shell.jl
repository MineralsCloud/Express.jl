module Shell

using Configurations: @option
using Mustache: render_from_file
using SimpleWorkflow: Script, ExternalAtomicJob, parallel

using ..Express: Action, Calculation

export ScriptTemplate, makescript

@option "template" struct ScriptTemplate
    file::String
    view::Dict
end

@option "script" struct ScriptConfig
    template::ScriptTemplate
end

function makescript(path, template::ScriptTemplate)
    str = render_from_file(template.file, template.view)
    if !isfile(path)
        if !isdir(dirname(path))
            mkpath(dirname(path))
        end
    else
        rm(path)
    end
    open(path, "w") do io
        write(io, str)
    end
    return Script(path)
end

function distprocs(nprocs, njobs)
    quotient, remainder = divrem(nprocs, njobs)
    if !iszero(remainder)
        @warn "The processes are not fully balanced! Consider the number of subjobs!"
    end
    return quotient
end

struct MakeCmd{T} <: Action{T} end

function buildjob(x::MakeCmd, inputs, args...; kwargs...)
    jobs = map(ExternalAtomicJob, x(inputs, 1; kwargs...))
    return parallel(jobs...)
end

end
