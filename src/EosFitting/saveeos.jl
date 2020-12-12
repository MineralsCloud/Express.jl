struct SaveEos{T} <: Action{T} end
function (::SaveEos{T})(path, eos::Parameters) where {T<:ScfOrOptim}
    ext = lowercase(extension(path))
    if ext == "jls"
        open(path, "w") do io
            serialize(io, eos)
        end
    elseif ext in ("json", "yaml", "yml", "toml")
        save(path, eos)
    else
        error("unsupported file extension `$ext`!")
    end
end
(x::SaveEos)(path, eos::EquationOfStateOfSolids) = x(path, getparam(eos))
