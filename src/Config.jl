module Config

using Configurations: @option
using Mustache: render_from_file
using SimpleWorkflow: Script

export ScriptTemplate, makescript

@option "script_template" struct ScriptTemplate
    file::String
    view::Dict
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

end
