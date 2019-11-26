module Cards

export card_helper

function card_helper end

module PWscf

using REPL.Terminals: TTYTerminal
using REPL.TerminalMenus: RadioMenu, request

using QuantumESPRESSO.Cards.PWscf: GammaPoint, MonkhorstPackGrid, KPointsCard

using ...Wizard: @c_str
using ..Cards

function Cards.card_helper(terminal::TTYTerminal, ::Type{T}) where {T<:KPointsCard}
    kpt_style = request(
        terminal,
        c"What k-point style do you want?"r,
        RadioMenu(["gamma", "automatic"]),
    )
    return if kpt_style == 1
        KPointsCard("gamma", GammaPoint())
    else  # "automatic"
        print(terminal, c"What 3-element k-point grid do you want (separated by spaces)?"r)
        grid = map(x -> parse(Int, x), split(readline(terminal), " ", keepempty = false))
        print(terminal, c"What 3-element k-point offsets do you want (separated by spaces)?"r)
        offsets = map(x -> parse(Int, x), split(readline(terminal), " ", keepempty = false))
        return KPointsCard("automatic", MonkhorstPackGrid(grid, offsets))
    end
end # function card_helper

end # module PWscf

end
