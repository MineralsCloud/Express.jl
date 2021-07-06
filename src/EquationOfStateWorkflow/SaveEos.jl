struct SaveEos{T} <: Action{T} end
function (::SaveEos{T})(file, eos::Parameters) where {T<:ScfOrOptim}
    ext = lowercase(extension(file))
    if ext == "jls"
        open(file, "w") do io
            serialize(io, eos)
        end
    elseif ext in ("json", "yaml", "yml", "toml")
        save(file, eos)
    else
        error("unsupported file extension `$ext`!")
    end
end
(x::SaveEos)(path, eos::EquationOfStateOfSolids) = x(path, getparam(eos))
