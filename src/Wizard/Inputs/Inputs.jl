module Inputs

using REPL
using REPL.Terminals
using REPL.TerminalMenus

using QuantumESPRESSO: to_qe
using QuantumESPRESSO.Namelists: Namelist
using QuantumESPRESSO.Namelists.PWscf
using QuantumESPRESSO.Inputs.PWscf: PWInput
using Rematch: @match
using Setfield: PropertyLens, set

using ..Namelists

export input_helper

function input_helper(terminal::TTYTerminal, ::Type{T}) where {T<:PWInput}
    control = namelist_helper(terminal, PWscf.ControlNamelist)
    system = namelist_helper(terminal, PWscf.SystemNamelist)
    electrons = namelist_helper(terminal, PWscf.ElectronsNamelist)
    ions = if control.calculation ∈ ("relax", "md", "vc-relax", "vc-md") 
        namelist_helper(terminal, PWscf.IonsNamelist)
    else
        PWscf.IonsNamelist()
    end
    cell = if control.calculation ∈ ("vc-relax", "vc-md")
        namelist_helper(terminal, PWscf.CellNamelist)
    else
        PWscf.CellNamelist()
    end
    # return T(
    #     control,
    #     system,
    #     electrons = ElectronsNamelist(),
    #     ions = IonsNamelist(),
    #     cell = CellNamelist(),
    # )
end # function input_helper

end
