module Namelists

using REPL
using REPL.Terminals
using REPL.TerminalMenus

using QuantumESPRESSO: to_qe
using QuantumESPRESSO.Namelists: Namelist
using QuantumESPRESSO.Namelists.PHonon
using QuantumESPRESSO.Namelists.PWscf
using Rematch: @match
using Setfield: PropertyLens, set

export namelist_helper

function setfield_helper(terminal::TTYTerminal, nml::T) where {T<:Namelist}
    while true
        print(terminal, to_qe(nml))
        isdone = pairs((false, true))[request(
            terminal,
            "We have generated an example `$(nameof(T))`. Want to change any field?",
            RadioMenu(["yes", "no"]),
        )]
        if !isdone
            print(terminal, "Type your field name: ")
            field = strip(readline(terminal)) |> Symbol
            if hasfield(T, field)
                print(terminal, "Type the value you want to set: ")
                nml = set(
                    nml,
                    PropertyLens{field}(),
                    parse(fieldtype(T, field), readline(terminal)),
                )
                continue
            end
            print(terminal, "Unknown field given! Try again!")
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
        "What exact calculation do you want to run?",
        RadioMenu(collect(values(calculations))),
    )]
    restart_mode =
        restart_modes[request(terminal, "Starting from scratch?", RadioMenu(["yes", "no"]))]
    print(terminal, "Convergence threshold on total energy (a.u): ")
    etot_conv_thr = parse(Float64, readline(terminal))
    print(terminal, "Convergence threshold on forces (a.u): ")
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
    print(terminal, "Please input the Bravais lattice index `ibrav`: ")
    ibrav = parse(Int, readline(terminal))
    print(terminal, "Please input a vector `celldm` 1-6: ")
    celldm = map(
        x -> parse(Float64, x),
        split(strip(readline(terminal), ['[', ']', '\n']), ",", keepempty = false),
    )
    print(terminal, "Please input the number of atoms in the unit cell `nat`: ")
    nat = parse(Int, readline(terminal))
    print(terminal, "Please input the number of types of atoms in the unit cell `ntyp`: ")
    ntyp = parse(Int, readline(terminal))
    print(
        terminal,
        "Please input the kinetic energy cutoff (Ry) for wavefunctions `ecutwfc`: ",
    )
    ecutwfc = parse(Float64, readline(terminal))
    print(
        terminal,
        "Please input the Kinetic energy cutoff (Ry) for charge density `ecutrho`: ",
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

function namelist_helper(terminal::TTYTerminal, ::Type{T}) where {T<:PHonon.PhNamelist}
    print(
        terminal, 
        "Please input the atomic mass [amu] of each atomic type `amass`: ",
    )
    amass = map(
        x -> parse(Float64, x),
        split(strip(readline(terminal), ['[', ']', '\n']), ",", keepempty = false),
    )
    epsil_pool = pairs((false, true))
    epsil = epsil_pool[request(
        terminal,
        "Please input the epsil `epsil`: ",
        RadioMenu(["false","true"]),
    )]
    q_in_band_form_pool = pairs((false, true))
    q_in_band_form = q_in_band_form_pool[request(
        terminal,
        "Please input the q_in_band_form `q_in_band_form`: ",
        RadioMenu(["false","true"]),
    )]
    print(terminal, "Please input parameters of the Monkhorst-Pack grid `nq` 1-3 in vector: ")
    nq1, nq2, nq3 = map(
        x -> parse(Float64, x),
        split(strip(readline(terminal), ['[', ']', '\n']), ",", keepempty = false),
    )
    print(terminal, "Please input parameters of the Monkhorst-Pack grid `nk` 1-3 in vector: ")
    nk1, nk2, nk3 = map(
        x -> parse(Float64, x),
        split(strip(readline(terminal), ['[', ']', '\n']), ",", keepempty = false),
    )
    print(terminal, "Please input offset `k` 1-3 in vector: ")
    k1, k2, k3 = map(
        x -> parse(Float64, x),
        split(strip(readline(terminal), ['[', ']', '\n']), ",", keepempty = false),
    )
    ph = T(
        amass = amass,
        epsil = epsil,
        q_in_band_form = q_in_band_form,
        nq1 = nq1,
        nq2 = nq2, 
        nq3 = nq3,
        nk1 = nk1, 
        nk2 = nk2, 
        nk3 = nk3,
        k1 = k1,
        k2 = k2,
        k3 = k3,
    )
    return setfield_helper(terminal, ph)
end # function namelist_helper
function namelist_helper(terminal::TTYTerminal, ::Type{T}) where {T<:PHonon.Q2rNamelist}
    print(
        terminal,
        "name of input dynamical matrices `fildyn`: ",
    )
    fildyn = readline(terminal)
    print(
        terminal,
        "name of output force constants `flfrc`: ",
    )
    flfrc = readline(terminal)
    zasr_pool = pairs(("no", "simple", "crystal", "one-dim", "zero-dim"))
    zasr = zasr_pool[request(
        terminal,
        "Please input the type of Acoustic Sum Rules used for the Born effective charges `zasr`: ",
        RadioMenu(collect(values(zasr_pool))),
    )]
    q2r = T(
        fildyn = fildyn,
        flfrc = flfrc,
        zasr = zasr,
    )
    return setfield_helper(terminal, q2r)
end # function namelist_helper
function namelist_helper(terminal::TTYTerminal, ::Type{T}) where {T<:PHonon.MatdynNamelist}
    dos_pool = pairs((false, true))
    dos = dos_pool[request(
        terminal,
        "Please input if calculate phonon Density of States `dos`: ",
        RadioMenu(["false","true"]),
    )]
    print(
        terminal,
        "Please input the energy step, in cm^(-1) `deltaE`: ",
    )
    deltaE = parse(Float64, readline(terminal))
    print(terminal, "Please input uniform q-point grid for DOS calculation `nk` 1-3 in vector: ")
    nk1, nk2, nk3 = map(
        x -> parse(Float64, x),
        split(strip(readline(terminal), ['[', ']', '\n']), ",", keepempty = false),
    )
    asr_pool = pairs(("no", "simple", "crystal", "one-dim", "zero-dim"))
    asr = asr_pool[request(
        terminal,
        "Please input the type of Acoustic Sum Rule `asr`: ",
        RadioMenu(collect(values(asr_pool))),
    )]
    print(
        terminal,
        "name of output force constants `flfrc`: ",
    )
    flfrc = readline(terminal)
    print(
        terminal,
        "name of input dynamical matrices `fildyn`: ",
    )
    fildyn = readline(terminal)
    print(
        terminal, 
        "Please input the masses of atoms in the supercell (a.m.u.) `amass`: ",
    )
    amass = map(
        x -> parse(Float64, x),
        split(strip(readline(terminal), ['[', ']', '\n']), ",", keepempty = false),
    )
    print(terminal, "Please input the number of atom types in the supercell `ntyp`: ")
    ntyp = parse(Int, readline(terminal))
    q_in_band_form_pool = pairs((false, true))
    q_in_band_form = q_in_band_form_pool[request(
        terminal,
        "Please input the q_in_band_form `q_in_band_form`: ",
        RadioMenu(["false","true"]),
    )]
    q_in_cryst_coord_pool = pairs((false, true))
    q_in_cryst_coord = q_in_cryst_coord_pool[request(
        terminal,
        "Please input the q_in_cryst_coord `q_in_cryst_coord`: ",
        RadioMenu(["false","true"]),
    )]
    nosym_pool = pairs((false, true))
    nosym = nosym_pool[request(
        terminal,
        "Please input if impose symmetry and time reversal `nosym`: ",
        RadioMenu(["false","true"]),
    )]
    matdyn = T(
        dos = dos,
        deltaE = deltaE,
        nk1 = nk1,
        nk2 = nk2,
        nk3 = nk3,
        asr = asr,
        flfrc = flfrc,
        fildyn = fildyn,
        amass = amass,
        ntyp = ntyp,
        q_in_band_form = q_in_band_form,
        q_in_cryst_coord = q_in_cryst_coord,
        nosym = nosym
    )
    return setfield_helper(terminal, matdyn)
end # function namelist_helper
function namelist_helper(terminal::TTYTerminal, ::Type{T}) where {T<:PHonon.DynmatNamelist}
    asr_pool = pairs(("no", "simple", "crystal", "one-dim", "zero-dim"))
    asr = asr_pool[request(
        terminal,
        "Please input the type of Acoustic Sum Rule `asr`: ",
        RadioMenu(collect(values(asr_pool))),
    )]
    print(
        terminal, 
        "Please input mass for each atom type `amass`: ",
    )
    amass = map(
        x -> parse(Float64, x),
        split(strip(readline(terminal), ['[', ']', '\n']), ",", keepempty = false),
    )
    dynmat = T(
        asr = asr,
        amass = amass,
    )
    return setfield_helper(terminal, dynmat)
end # function namelist_helper
end
