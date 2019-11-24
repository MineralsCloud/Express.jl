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

function pwscf_helper(terminal::TTYTerminal)
    calculations = pairs(("scf", "nscf", "bands", "relax", "md", "vc-relax", "vc-md"))
    restart_modes = pairs(("from_scratch", "restart"))
    calculation = calculations[request(terminal,
        "What exact calculation do you want to run?",
        RadioMenu(collect(values(calculations))),
    )]
    restart_mode = restart_modes[request(terminal,
        "Starting from scratch?",
        RadioMenu(["yes", "no"]),
    )]
    print(terminal, "Convergence threshold on total energy (a.u): ")
    etot_conv_thr = parse(Float64, readline(terminal))
    print(terminal, "Convergence threshold on forces (a.u): ")
    forc_conv_thr = parse(Float64, readline(terminal))
    control = PWscf.ControlNamelist(
        calculation = calculation,
        restart_mode = restart_mode,
        etot_conv_thr = etot_conv_thr,
        forc_conv_thr = forc_conv_thr,
    )
    while true
        print(terminal, to_qe(control))
        isdone = pairs((false, true))[request(terminal,
            "We have generated an example `ControlNamelist`. Want to change any field?",
            RadioMenu(["yes", "no"]),
        )]
        if !isdone
            print(terminal, "Type your field name: ")
            field = strip(readline(terminal)) |> Symbol
            if hasfield(PWscf.ControlNamelist, field)
                print(terminal, "Type the value you want to set: ")
                control = set(control, PropertyLens{field}(), parse(fieldtype(PWscf.ControlNamelist, field), readline(terminal)))
                continue
            end
            print(terminal, "Unknown field given! Try again!")
        end
        break
    end
    return control
end # function pwscf_helper

end
