module Config

using Configurations: @option
using Mustache: render_from_file
using SimpleWorkflow: Script

export ShellTemplate, makescript

@option "shell_template" struct ShellTemplate
    file::String
    view::Dict
end

function makescript(path, template::ShellTemplate)
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
