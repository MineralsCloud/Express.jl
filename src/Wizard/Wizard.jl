module Wizard

using REPL
using REPL.Terminals
using REPL.TerminalMenus

using Compat: isnothing
using Crayons.Box: GREEN_FG
using EquationsOfState.Collections
using JLD2: jldopen
using Parameters: @with_kw
using QuantumESPRESSO: to_qe
using QuantumESPRESSO.Namelists.PWscf
using QuantumESPRESSO.Inputs.PHonon: PhInput, Q2rInput, MatdynInput, DynmatInput
using QuantumESPRESSO.Inputs.PWscf: PWInput
using QuantumESPRESSO.Inputs.CP: CPInput
using QuantumESPRESSOHelpers.Inputs: input_builder
using Setfield: PropertyLens, set

using ..EquationOfStateFitting: update_alat_press

export run_wizard

abstract type QuantumESPRESSOCalculation end
struct PWscfCalculation <: QuantumESPRESSOCalculation end
struct PHononCalculation <: QuantumESPRESSOCalculation end
struct CPCalculation <: QuantumESPRESSOCalculation end

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
            printstyled(state.out, msg, bold = true, color = :red)
        else
            bt = catch_backtrace()
            Base.showerror(stderr, err, bt)
            println(state.out, "\n")
        end
        return state
    end

    # We did it!
    save_last_wizard_state(state)

    println(state.out, GREEN_FG("Wizard Complete. Press any key to exit..."))
    read(state.in, Char)

    return state
end # function run_wizard

input_type(::PWscfCalculation) = PWInput
input_type(::PHononCalculation) = PhInput
input_type(::CPCalculation) = CPInput

step(i::Integer, state::WizardState) = step(Val(i), state)
function step(::Val{1}, state::WizardState)
    terminal = TTYTerminal("xterm", state.in, state.out, state.out)
    calculation = pairs((PWscfCalculation(), PHononCalculation(), CPCalculation()))[request(
        terminal,
        GREEN_FG("What calculation do you want to run?") |> string,
        RadioMenu(["scf", "phonon", "CPMD"]),
    )]
    state.result = input_builder(terminal, input_type(calculation))
    return state
end # function step
function step(::Val{2}, state::WizardState)
    return step(Val(2), state, PWscfCalculation())
end # function step
function step(::Val{2}, state::WizardState, ::PWscfCalculation)
    terminal = TTYTerminal("xterm", state.in, state.out, state.out)
    eos_symb = Symbol(request(
        terminal,
        GREEN_FG("What EOS do you want to use for fitting?") |> string,
        RadioMenu([
            "Murnaghan",
            "BirchMurnaghan2nd",
            "BirchMurnaghan3rd",
            "BirchMurnaghan4th",
            "PoirierTarantola2nd",
            "PoirierTarantola3rd",
            "PoirierTarantola4th",
            "Vinet",
        ]),
    ))
    println(terminal, GREEN_FG("Please input parameters for this EOS (separated by spaces):") |> string)
    eos = eval(eos_symb)(map(
        x -> parse(Float64, x),
        split(readline(terminal), " ", keepempty = false),
    ))
    println(terminal, GREEN_FG("Please input pressures you want to test on (separated by spaces):") |> string)
    pressures =
        map(x -> parse(Float64, x), split(readline(terminal), " ", keepempty = false))
    inputs_from_template =
        [update_alat_press(state.result, eos, pressure) for pressure in pressures]

end # function step

end
