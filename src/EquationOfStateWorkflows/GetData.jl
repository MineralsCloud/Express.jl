struct GetData{T} <: Action{T} end
function (::GetData{T})(outputs) where {T<:ScfOrOptim}
    raw = (parseoutput(T())(output) for output in outputs)  # `ntuple` cannot work with generators
    return collect(Iterators.filter(!isnothing, raw))  # A vector of pairs
end
function (x::GetData{T})(file, outputs) where {T}
    data = x(outputs)
    dict = Dict(
        "volume" => (ustrip ∘ first).(data),
        "energy" => (ustrip ∘ last).(data),
        "vunit" => string(unit(first(data).first)),
        "eunit" => string(unit(first(data).second)),
    )
    ext = lowercase(extension(file))
    if ext == "jls"
        open(file, "w") do io
            serialize(io, dict)
        end
    elseif ext in ("json", "yaml", "yml", "toml")
        save(file, dict)
    else
        error("unsupported file extension `$ext`!")
    end
end

function parseoutput end
