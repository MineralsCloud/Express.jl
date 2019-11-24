module Wizard

using REPL
using REPL.Terminals
using REPL.TerminalMenus

using QuantumESPRESSO: to_qe
using QuantumESPRESSO.Namelists.PWscf
using Rematch: @match
using Setfield: PropertyLens, set

export run_wizard

Base.@kwdef mutable struct WizardState
    step::Int = 1
    ins::IO = stdin
    outs::IO = stdout
end

function run_wizard(state::Union{Nothing,WizardState} = nothing)
    terminal = TTYTerminal("xterm", state.ins, state.outs, state.outs)
    choice = request(terminal,
        "Would you like to resume the previous incomplete wizard run?",
        RadioMenu([
            "Resume previous run",
            "Start from scratch",
        ]),
    )
    # if choice == 1
    #     return state
    # else
    #     return WizardState()
    # end
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
end # function run_wizard

function pwscf_helper(terminal)
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
