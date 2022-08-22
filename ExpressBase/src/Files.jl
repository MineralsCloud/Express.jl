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

format(::Val{:json}) = DataFormat{:JSON}()
format(::Val{:yaml}) = DataFormat{:YAML}()
format(::Val{:yml}) = DataFormat{:YAML}()
format(::Val{:toml}) = DataFormat{:TOML}()
@valsplit format(Val(ext::Symbol)) = throw(UnsupportedExtensionError(string(ext)))

# See https://github.com/JuliaIO/FileIO.jl/blob/b779539/src/types.jl#L16-L18
macro format_str(s)
    return :(DataFormat{$(Expr(:quote, Symbol(s)))})
end

"""
    save(file, data)

Save `data` to `file`.

By now, `YAML`, `JSON`, and `TOML` formats are supported. The format is recognized by `file` extension.

If `data` is a `Dict`, its keys should be `String`s so that `load` can return the same `data`.

!!! warning
    Allowed `data` types can be referenced in [`JSON.jl` documentation](https://github.com/JuliaIO/JSON.jl/blob/master/README.md)
    and [`YAML.jl` documentation](https://github.com/JuliaData/YAML.jl/blob/master/README.md).
    For `TOML` format, only `AbstractDict` type is allowed.
"""
function save(file, data)
    path, ext = expanduser(file), extension(file)
    save(File{format(ext)}(path), data)
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

By now, `YAML`, `JSON`, and `TOML` formats are supported. The format is recognized by `file` extension.
"""
function load(file)
    path, ext = expanduser(file), extension(file)
    return load(File{format(ext)}(path))
end
load(path::File{format"JSON"}) = JSON.parsefile(path)
function load(path::File{format"TOML"})
    open(path, "r") do io
        return TOML.parse(io)
    end
end
function load(path::File{format"YAML"})
    open(path, "r") do io
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
Get the extension from `file`. Return an empty string if no extension is found.
"""
function extension(file)
    ext = splitext(file)[2]  # `splitext` works with `FilePathsBase.AbstractPath` since version `v0.7.0`.
    return isempty(ext) ? "" : lowercase(ext[2:end])
end

function Base.show(io::IO, error::UnsupportedExtensionError)
    return print(io, "unsupported extension `.", error.ext, "`!")
end

end
