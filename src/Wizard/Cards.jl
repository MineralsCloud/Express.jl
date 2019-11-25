module Cards

using REPL
using REPL.Terminals
using REPL.TerminalMenus

using QuantumESPRESSO: to_qe
using QuantumESPRESSO.Cards.PWscf: GammaPoint, MonkhorstPackGrid, KPointsCard
using Rematch: @match
using Setfield: PropertyLens, set

using ..Wizard: @c_str

function card_helper(terminal::TTYTerminal, ::Type{T}) where {T<:PWscf.KPointsCard}
    kpt_style = request(
        terminal,
        c"What k-point style do you want?"r,
        RadioMenu(["gamma", "automatic"]),
    )
    return if kpt_style == 1
        KPointsCard("gamma", GammaPoint())
    else  # "automatic"
        print(terminal, "What 3-element k-point grid do you want (separated by spaces)?")
        grid = map(x -> parse(Int, x), split(readline(terminal), " ", keepempty = false))
        print(terminal, "What 3-element k-point offsets do you want (separated by spaces)?")
        offsets = map(x -> parse(Int, x), split(readline(terminal), " ", keepempty = false))
        return KPointsCard("automatic", MonkhorstPackGrid(grid, offsets))
    end
end # function card_helper

end
