module Inputs

using REPL
using REPL.Terminals
using REPL.TerminalMenus

using QuantumESPRESSO: to_qe
using QuantumESPRESSO.Namelists: Namelist
using QuantumESPRESSO.Namelists.PWscf
using Rematch: @match
using Setfield: PropertyLens, set

using ..Namelists

export pwscf_helper

function pwscf_helper(terminal::TTYTerminal)
    control = namelist_helper(terminal, PWscf.ControlNamelist)
    system = namelist_helper(terminal, PWscf.SystemNamelist)
    return (control, system)
end # function pwscf_helper

end
