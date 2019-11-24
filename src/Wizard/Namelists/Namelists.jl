module Namelists

using REPL
using REPL.Terminals
using REPL.TerminalMenus

using QuantumESPRESSO: to_qe
using QuantumESPRESSO.Namelists: Namelist
using QuantumESPRESSO.Namelists.PWscf
using Rematch: @match
using Setfield: PropertyLens, set

using ..Wizard: @c_str

export namelist_helper

function setfield_helper(terminal::TTYTerminal, nml::T) where {T<:Namelist}
    while true
        print(terminal, to_qe(nml))
        isdone = pairs((false, true))[request(
            terminal,
            c"We have generated an example `$(nameof(T))`. Want to change any field?"r,
            RadioMenu(["yes", "no"]),
        )]
        if !isdone
            print(terminal, c"Type your field name: "r)
            field = strip(readline(terminal)) |> Symbol
            if hasfield(T, field)
                print(terminal, c"Type the value you want to set: "r)
                nml = set(
                    nml,
                    PropertyLens{field}(),
                    parse(fieldtype(T, field), readline(terminal)),
                )
                continue
            end
            print(terminal, c"Unknown field given! Try again!"r)
        end
        break
    end
    return nml
end # function change_field_helper

function namelist_helper(terminal::TTYTerminal, ::Type{T}) where {T<:PWscf.ControlNamelist}
    calculations = pairs(("scf", "nscf", "bands", "relax", "md", "vc-relax", "vc-md"))
    restart_modes = pairs(("from_scratch", "restart"))
    calculation = calculations[request(
        terminal,
        c"What exact calculation do you want to run?"r,
        RadioMenu(collect(values(calculations))),
    )]
    restart_mode = restart_modes[request(
        terminal,
        c"Starting from scratch?"r,
        RadioMenu(["yes", "no"]),
    )]
    print(terminal, c"Convergence threshold on total energy (a.u): "r)
    etot_conv_thr = parse(Float64, readline(terminal))
    print(terminal, c"Convergence threshold on forces (a.u): "r)
    forc_conv_thr = parse(Float64, readline(terminal))
    control = T(
        calculation = calculation,
        restart_mode = restart_mode,
        etot_conv_thr = etot_conv_thr,
        forc_conv_thr = forc_conv_thr,
    )
    return setfield_helper(terminal, control)
end # function namelist_helper
function namelist_helper(terminal::TTYTerminal, ::Type{T}) where {T<:PWscf.SystemNamelist}
    print(terminal, c"Please input the Bravais lattice index `ibrav`: "r)
    ibrav = parse(Int, readline(terminal))
    print(terminal, c"Please input a vector `celldm` 1-6: "r)
    celldm = map(
        x -> parse(Float64, x),
        split(strip(readline(terminal), ['[', ']', '\n']), ",", keepempty = false),
    )
    print(terminal, c"Please input the number of atoms in the unit cell `nat`: "r)
    nat = parse(Int, readline(terminal))
    print(terminal, c"Please input the number of types of atoms in the unit cell `ntyp`: "r)
    ntyp = parse(Int, readline(terminal))
    print(
        terminal,
        c"Please input the kinetic energy cutoff (Ry) for wavefunctions `ecutwfc`: "r,
    )
    ecutwfc = parse(Float64, readline(terminal))
    print(
        terminal,
        c"Please input the Kinetic energy cutoff (Ry) for charge density `ecutrho`: "r,
    )
    ecutrho = parse(Float64, readline(terminal))
    system = T(
        ibrav = ibrav,
        celldm = celldm,
        nat = nat,
        ntyp = ntyp,
        ecutwfc = ecutwfc,
        ecutrho = ecutrho,
    )
    return setfield_helper(terminal, system)
end # function namelist_helper
function namelist_helper(
    terminal::TTYTerminal,
    ::Type{T},
) where {T<:PWscf.ElectronsNamelist}
    print(
        terminal,
        "Please input the convergence threshold for selfconsistency `conv_thr`: ",
    )
    conv_thr = parse(Float64, readline(terminal))
    diagonalizations = pairs(("david", "cg", "cg-serial", "david-serial"))
    diagonalization = diagonalizations[request(
        terminal,
        "Please input the diagonalization method `diagonalization`: ",
        RadioMenu(collect(values(diagonalizations))),
    )]
    electrons = T(conv_thr = conv_thr, diagonalization = diagonalization)
    return setfield_helper(terminal, electrons)
end # function namelist_helper
function namelist_helper(terminal::TTYTerminal, ::Type{T}) where {T<:PWscf.IonsNamelist}
    ion_dynamics_pool =
        pairs(("none", "bfgs", "damp", "verlet", "langevin", "langevin-smc", "beeman"))
    ion_dynamics = ion_dynamics_pool[request(
        terminal,
        "Please input the type of ionic dynamics `ion_dynamics`: ",
        RadioMenu(collect(values(ion_dynamics_pool))),
    )]
    ion_temperature_pool = (
        "rescaling",
        "rescale-v",
        "rescale-T",
        "reduce-T",
        "berendsen",
        "andersen",
        "initial",
        "not_controlled",
    )
    ion_temperature = ion_temperature_pool[request(
        terminal,
        "Please input the ions temperature `ion_temperature`: ",
        RadioMenu(collect(values(ion_temperature_pool))),
    )]
    ions = T(ion_dynamics = ion_dynamics, ion_temperature = ion_temperature)
    return setfield_helper(terminal, ions)
end # function namelist_helper
function namelist_helper(terminal::TTYTerminal, ::Type{T}) where {T<:PWscf.CellNamelist}
    cell_dynamics_pool = pairs(("none", "sd", "damp-pr", "damp-w", "bfgs", "pr", "w"))
    cell_dynamics = cell_dynamics_pool[request(
        terminal,
        "Please input the type of dynamics for the cell `cell_dynamics`: ",
        RadioMenu(collect(values(cell_dynamics_pool))),
    )]
    print(
        terminal,
        "Please input the target pressure [KBar] in a variable-cell md or relaxation run `press`: ",
    )
    press = parse(Float64, readline(terminal))
    print(
        terminal,
        "Please input the fictitious cell mass [amu] for variable-cell simulations `wmass`: ",
    )
    wmass = parse(Float64, readline(terminal))
    print(
        terminal,
        "Please input the Convergence threshold on the pressure for variable cell `press_conv_thr`: ",
    )
    press_conv_thr = parse(Float64, readline(terminal))
    cell = T(
        cell_dynamics = cell_dynamics,
        press = press,
        wmass = wmass,
        press_conv_thr = press_conv_thr,
    )
    return setfield_helper(terminal, cell)
end # function namelist_helper
end
