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

# Referenced from https://github.com/JuliaPackaging/BinaryBuilder.jl/blob/0eece73/src/wizard/state.jl
function save_last_wizard_state(state::WizardState)
    dir = joinpath(dirname(@__DIR__), "wizard_state")
    isdir(dir) || mkdir(dir)
    jldopen(joinpath(dir, "wizard.state"), "w") do f
        serialize(f, state)
    end
    return state
end

# Referenced from https://github.com/JuliaPackaging/BinaryBuilder.jl/blob/0eece73/src/wizard/state.jl
function load_last_wizard_state()
    wizard_state_dir = joinpath(dirname(@__DIR__), "wizard_state")

    # If no state dir exists, early-exit
    if isnothing(wizard_state_dir)
        return WizardState()
    end

    try
        state = jldopen(joinpath(wizard_state_dir, "wizard.state"), "r") do f
            return deserialize(f)
        end

        # Looks like we had an incomplete build; ask the user if they want to continue
        if !(state.step ∈ (0, 1))  # 0: end, 1: start
            terminal = TTYTerminal("xterm", state.ins, state.outs, state.outs)
            choice = request(terminal,
                "Would you like to resume the previous incomplete wizard run?",
                RadioMenu([
                    "Resume previous run",
                    "Start from scratch",
                ]),
            )

            if choice == 1
                return state
            else
                return WizardState()
            end
        end
    catch e
        if isa(e, InterruptException)
            rethrow(e)
        end
    end

    # Either something went wrong, or there was nothing interesting stored.
    # Either way, just return a blank slate.
    return WizardState()
end
