using JSON
using YAML

function save(filepath, data)
    ext = extension(filepath)
    if ext ∈ ("yaml", "yml")
        YAML.write_file(expanduser(filepath), data)
    elseif ext == "json"
        open(expanduser(filepath), "w") do io
            JSON.print(io, data)
        end
    else
        error("unknown file extension `$ext`!")
    end
end # function save

function load(filepath)
    ext = extension(filepath)
    if ext ∈ ("yaml", "yml")
        return open(expanduser(filepath), "r") do io
            YAML.load(io)
        end
    elseif ext == "json"
        return JSON.parsefile(expanduser(filepath))
    else
        error("unknown file extension `$ext`!")
    end
end # function load

function extension(filepath)  # From https://github.com/rofinn/FilePathsBase.jl/blob/af850a4/src/path.jl#L331-L340
    name = basename(filepath)
    tokenized = split(name, '.')
    if length(tokenized) > 1
        return lowercase(tokenized[end])
    else
        return ""
    end
end
