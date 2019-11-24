module Wizard

using REPL
using REPL.Terminals
using REPL.TerminalMenus

using Compat: isnothing
using JLD2: jldopen
using Parameters: @with_kw
using QuantumESPRESSO: to_qe
using QuantumESPRESSO.Namelists.PWscf
using QuantumESPRESSO.Inputs.PHonon: PhInput, Q2rInput, MatdynInput, DynmatInput
using QuantumESPRESSO.Inputs.PWscf: PWInput
using Rematch: @match
using Setfield: PropertyLens, set

export run_wizard

include("utils.jl")
include("state.jl")
include("Namelists/Namelists.jl")
include("Inputs/Inputs.jl")
using .Inputs: input_helper

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
            printstyled(state.outs, msg, bold = true, color = :red)
        else
            bt = catch_backtrace()
            Base.showerror(stderr, err, bt)
            println(state.outs, "\n")
        end
        return state
    end

    # We did it!
    save_last_wizard_state(state)

    println(state.outs, c"Wizard Complete. Press any key to exit..."g)
    read(state.ins, Char)

    return state
end # function run_wizard

step(i::Integer, state::WizardState) = step(Val(i), state)
function step(::Val{1}, state::WizardState)
    terminal = TTYTerminal("xterm", state.ins, state.outs, state.outs)
    choice = request(
        terminal,
        c"What calculation do you want to run?"r,
        RadioMenu(["scf", "phonon", "CPMD"]),
    )
    @match choice begin
        1 => input_helper(terminal, PWInput)
        2 => input_helper(terminal, PhInput)
    end
end # function step
function step(::Val{2}, state::WizardState)

end # function step

end
