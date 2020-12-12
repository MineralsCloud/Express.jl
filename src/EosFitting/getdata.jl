struct GetData{T} <: Action{T} end
function (::GetData{T})(outputs) where {T<:ScfOrOptim}
    raw = (load(parseoutput(T()), output) for output in outputs)  # `ntuple` cannot work with generators
    return collect(Iterators.filter(!isnothing, raw))  # A vector of pairs
end

function parseoutput end
