module Shell

using Configurations: @option
using Mustache: render_from_file
using SimpleWorkflow: InternalAtomicJob, Script, parallel

export ScriptTemplate, makescript, @intjob

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

macro intjob(ex)  # See https://github.com/JuliaLang/julia/blob/ab5853f/base/task.jl#L111-L113
    return :(InternalAtomicJob(() -> $(esc(ex))))
end

end
