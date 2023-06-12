module Files

using JSON: JSON
using TOML: TOML
using ValSplit: @valsplit
using YAML: YAML

export load, save, of_format, @format_str

# See https://github.com/JuliaIO/FileIO.jl/blob/b779539/src/types.jl#L14
struct DataFormat{T} end

# See https://github.com/JuliaIO/FileIO.jl/blob/b779539/src/types.jl#L38-L44
struct File{D<:DataFormat,N}
    filename::N
end
File{F}(file::File{F}) where {F<:DataFormat} = file
File{F}(file::AbstractString) where {F<:DataFormat} = File{F,String}(String(file))
File{F}(file) where {F<:DataFormat} = File{F,typeof(file)}(file)
function File(file)
    ext = extension(file)
    fmt = format(Symbol(ext))
    return File{fmt}(file)
end

struct UnsupportedExtensionError <: Exception
    ext::String
end

# See https://github.com/JuliaIO/FileIO.jl/blob/b779539/src/types.jl#L16-L18
macro format_str(s)
    return :(DataFormat{$(Expr(:quote, Symbol(s)))})
end

format(::Val{:json}) = format"JSON"
format(::Val{:yaml}) = format"YAML"
format(::Val{:yml}) = format"YAML"
format(::Val{:toml}) = format"TOML"
@valsplit format(Val(ext::Symbol)) = throw(UnsupportedExtensionError(string(ext)))

"""
    save(file, data)

Save `data` to `file`.

By now, `YAML`, `JSON`, and `TOML` formats are supported. The format is recognized by the file extension.

If `data` is a `Dict`, its keys should be `String`s so that `load` can return the same `data`.

!!! warning
    Allowed `data` types can be referenced in [`JSON.jl` documentation](https://github.com/JuliaIO/JSON.jl/blob/master/README.md)
    and [`YAML.jl` documentation](https://github.com/JuliaData/YAML.jl/blob/master/README.md).
    For `TOML` format, only `AbstractDict` type is allowed.
"""
function save(file, data)
    save(File(expanduser(file)), data)
    return nothing
end
function save(file::File{format"JSON"}, data)
    open(file, "w") do io
        JSON.print(io, data)
    end
end
function save(file::File{format"TOML"}, data)
    open(file, "w") do io
        TOML.print(io, data)
    end
end
function save(file::File{format"YAML"}, data)
    open(file, "w") do io
        YAML.write(io, data, "")
    end
end

"""
    load(file)

Load data from `file` to a `Dict`.

By now, `YAML`, `JSON`, and `TOML` formats are supported. The format is recognized by the file extension.
"""
load(file) = load(File(expanduser(file)))
function load(file::File{format"JSON"})
    open(file, "r") do io
        str = read(io, String)
        return JSON.parse(str)
    end
end
function load(file::File{format"TOML"})
    open(file, "r") do io
        return TOML.parse(io)
    end
end
function load(file::File{format"YAML"})
    open(file, "r") do io
        dict = YAML.load(io)
        return JSON.parse(JSON.json(dict))  # To keep up with JSON & TOML results
    end
end

"""
    of_format(destination, source)

Convert `source` to the format of `destination`. Similar to `oftype`.
"""
function of_format(destination, source)
    data = load(source)
    return save(destination, data)
end

"""
    extension(file)

Get the extension from `file`. Return an empty `String` if no extension is found.
"""
function extension(file)
    ext = splitext(file)[2]  # `splitext` works with `FilePathsBase.AbstractPath` since version `v0.7.0`.
    return isempty(ext) ? "" : lowercase(ext[2:end])
end

# See https://github.com/JuliaIO/FileIO.jl/blob/b779539/src/types.jl#L121-L124
Base.open(@nospecialize(file::File), @nospecialize(args...)) = open(file.filename, args...)

Base.close(@nospecialize(file::File)) = close(file.filename)

Base.read(@nospecialize(file::File), @nospecialize(args...)) = read(file.filename, args...)

function Base.show(io::IO, error::UnsupportedExtensionError)
    return print(io, "unsupported extension `.", error.ext, "`!")
end

end
