module Config

export ConfigFile

struct ConfigFile{T}
    source::T
end

Base.open(file::ConfigFile, args...; kwargs...) = open(file.source, args...; kwargs...)
Base.open(f::Function, file::ConfigFile, args...; kwargs...) =
    open(f, file.source, args...; kwargs...)

include("sp.jl")
include("io.jl")

end
