module Inputs

export input_helper

function input_helper end

module PWscf

using REPL.Terminals: TTYTerminal
using REPL.TerminalMenus: RadioMenu, request

using QuantumESPRESSO.Namelists.PWscf:
    ControlNamelist, SystemNamelist, ElectronsNamelist, IonsNamelist, CellNamelist
using QuantumESPRESSO.Cards.PWscf: AtomicSpecies, AtomicSpeciesCard, AtomicPosition, AtomicPositionsCard, KPointsCard
using QuantumESPRESSO.Inputs.PWscf: PWInput

using ...Namelists: namelist_helper
using ...Cards: card_helper
using ..Inputs

function Inputs.input_helper(terminal::TTYTerminal, ::Type{T}) where {T<:PWInput}
    control = namelist_helper(terminal, ControlNamelist)
    system = namelist_helper(terminal, SystemNamelist)
    electrons = namelist_helper(terminal, ElectronsNamelist)
    ions = if control.calculation ∈ ("relax", "md", "vc-relax", "vc-md")
        namelist_helper(terminal, IonsNamelist)
    else
        IonsNamelist()
    end
    cell = if control.calculation ∈ ("vc-relax", "vc-md")
        namelist_helper(terminal, CellNamelist)
    else
        CellNamelist()
    end
    k_points = card_helper(terminal, KPointsCard)
    return T(
        control = control,
        system = system,
        electrons = electrons,
        ions = ions,
        cell = cell,
        atomic_species = AtomicSpeciesCard(AtomicSpecies[]),
        atomic_positions = AtomicPositionsCard("alat", AtomicPosition[]),
        k_points = k_points,
        cell_parameters = nothing,
    )
end # function input_helper

end # module PWscf

module PHonon

using REPL.Terminals: TTYTerminal
using REPL.TerminalMenus: RadioMenu, request

using QuantumESPRESSO.Namelists.PHonon:
    PhNamelist, Q2rNamelist, MatdynNamelist, DynmatNamelist
using QuantumESPRESSO.Inputs.PHonon: PhInput, Q2rInput, MatdynInput, DynmatInput

using ...Namelists: namelist_helper
using ..Inputs

function Inputs.input_helper(terminal::TTYTerminal, ::Type{T}) where {T<:PhInput}
    ph = namelist_helper(terminal, PhNamelist)
end
function Inputs.input_helper(terminal::TTYTerminal, ::Type{T}) where {T<:Q2rInput}
    return Q2rInput(namelist_helper(terminal, Q2rNamelist))
end
function Inputs.input_helper(terminal::TTYTerminal, ::Type{T}) where {T<:MatdynInput}
    ph = namelist_helper(terminal, MatdynNamelist)
end
function Inputs.input_helper(terminal::TTYTerminal, ::Type{T}) where {T<:DynmatInput}
    return DynmatInput(namelist_helper(terminal, DynmatNamelist))
end

end # module PHonon

end
