module Wizard

using REPL
using REPL.Terminals
using REPL.TerminalMenus

using Compat: isnothing
using JLD2: jldopen
using Parameters: @with_kw
using QuantumESPRESSO: to_qe
using QuantumESPRESSO.Namelists.PWscf
using Rematch: @match
using Setfield: PropertyLens, set

export run_wizard

include("state.jl")
include("Namelists/Namelists.jl")
include("Inputs/Inputs.jl")
using .Inputs: pwscf_helper

# Referenced from https://github.com/JuliaPackaging/BinaryBuilder.jl/blob/0eece73/src/wizard/state.jl
function run_wizard(state::Union{Nothing,WizardState} = nothing)
    if isnothing(state)
        # If we weren't given a state, check to see if we'd like to resume a
        # previous run or start from scratch again.
        state = load_last_wizard_state()
    end

    try
        while state.step != 0
            if state.step == 1
                step(1, state)
                state.step = 2
            elseif state.step == 2
                step(2, state)
                state.step = 0  # End step
            end

            # Save it every step along the way
            save_last_wizard_state(state)
        end
    catch err
        # If anything goes wrong, immediately save the current wizard state
        save_last_wizard_state(state)
        if isa(err, InterruptException)
            msg = "\n\nWizard stopped, use run_wizard() to resume.\n\n"
            printstyled(state.outs, msg, bold=true, color=:red)
        else
            bt = catch_backtrace()
            Base.showerror(stderr, err, bt)
            println(state.outs, "\n")
        end
        return state
    end

    # We did it!
    save_last_wizard_state(state)

    println(state.outs, "\nWizard Complete. Press any key to exit...")
    read(state.ins, Char)

    return state
end # function run_wizard

step(i::Integer, state::WizardState) = step(Val(i), state)
function step(::Val{1}, state::WizardState)
    terminal = TTYTerminal("xterm", state.ins, state.outs, state.outs)
    choice = request(terminal,
        "What calculation do you want to run?",
        RadioMenu([
            "scf",
            "phonon",
            "CPMD",
        ]),
    )
    @match choice begin
        1 => pwscf_helper(terminal)
    end
end # function step
function step(::Val{2}, state::WizardState)
    
end # function step

end
