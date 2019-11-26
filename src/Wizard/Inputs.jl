module Inputs

export input_helper

function input_helper end

module PWscf

using REPL.Terminals: TTYTerminal
using REPL.TerminalMenus: RadioMenu, request

using QuantumESPRESSO: asfieldname
using QuantumESPRESSO.Namelists.PWscf:
    ControlNamelist, SystemNamelist, ElectronsNamelist, IonsNamelist, CellNamelist
using QuantumESPRESSO.Cards.PWscf: AtomicSpecies, AtomicSpeciesCard, AtomicPosition, AtomicPositionsCard, KPointsCard, CellParametersCard
using QuantumESPRESSO.Inputs.PWscf: PWInput

using ...Namelists: namelist_helper
using ...Cards: card_helper
using ...Wizard: @c_str
using ..Inputs

function Inputs.input_helper(terminal::TTYTerminal, ::Type{T}) where {T<:PWInput}
    fields = Dict{Symbol,Any}()
    for S in (ControlNamelist, SystemNamelist, ElectronsNamelist)
        haserror = true
        while haserror
            try
                push!(fields, asfieldname(S) => namelist_helper(terminal, S))
                haserror = false
            catch e
                isa(e, InterruptException) && rethrow(e)
                println(terminal, c"Something wrong happens, try again!"g)
            end
        end
    end
    if fields[:control].calculation ∈ ("relax", "md", "vc-relax", "vc-md")
        haserror = true
        while haserror
            try
                push!(fields, asfieldname(IonsNamelist) => namelist_helper(terminal, IonsNamelist))
                haserror = false
            catch e
                isa(e, InterruptException) && rethrow(e)
                println(terminal, c"Something wrong happens, try again!"g)
            end
        end
    else
        IonsNamelist()
    end
    if fields[:control].calculation ∈ ("vc-relax", "vc-md")
        haserror = true
        while haserror
            try
                push!(fields, asfieldname(CellNamelist) => namelist_helper(terminal, CellNamelist))
                haserror = false
            catch e
                isa(e, InterruptException) && rethrow(e)
                println(terminal, c"Something wrong happens, try again!"g)
            end
        end
    else
        CellNamelist()
    end
    haserror = true
    while haserror
        try
            push!(fields, asfieldname(KPointsCard) => card_helper(terminal, KPointsCard))
            haserror = false
        catch e
            isa(e, InterruptException) && rethrow(e)
            println(terminal, c"Something wrong happens, try again!"g)
        end
    end
    push!(fields, asfieldname(AtomicSpeciesCard) => AtomicSpeciesCard(AtomicSpecies[]))
    push!(fields, asfieldname(AtomicPositionsCard) => AtomicPositionsCard("alat", AtomicPosition[]))
    push!(fields, asfieldname(CellParametersCard) => nothing)
    return T(; fields...)
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
