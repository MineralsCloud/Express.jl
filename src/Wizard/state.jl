# Referenced from https://github.com/JuliaPackaging/BinaryBuilder.jl/blob/0eece73/src/wizard/state.jl
Base.@kwdef mutable struct WizardState
    step::Int = 1
    ins::IO = stdin
    outs::IO = stdout
end

# Serialize a WizardState out into a JLD2 dictionary-like object
function serialize(io, x::WizardState)
    for field in fieldnames(typeof(x))
        io[string(field)] = getproperty(x, field)
    end
    # For non-serializable fields (such as `x.ins` and `x.outs`) we just recreate them in unserialize().
end

function deserialize(io)
    x = WizardState()
    for field in fieldnames(typeof(x))
        setproperty!(x, field, io[string(field)])
    end
    # Manually recreate `ins` and `outs`.  Note that this just sets them to their default values
    x.ins = stdin
    x.outs = stdout
    return x
end
