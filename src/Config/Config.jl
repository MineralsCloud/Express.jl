module Config

export ConfigFile

struct ConfigFile{T}
    source::T
end

Base.open(file::ConfigFile, args...; kwargs...) = open(file.source, args...; kwargs...)
function Base.open(f::Function, file::ConfigFile, args...; kwargs...)
    return open(f, file.source, args...; kwargs...)
end

include("sp.jl")
include("io.jl")

end
