module Namelists

using REPL.Terminals: TTYTerminal
using REPL.TerminalMenus: RadioMenu, request

using QuantumESPRESSO: to_qe
using QuantumESPRESSO.Namelists: Namelist
using Setfield: PropertyLens, set

using ..Wizard: color_string, @c_str

export namelist_helper

function namelist_helper end

# This is a helper function and should not be exported.
function setfield_helper(terminal::TTYTerminal, nml::T) where {T<:Namelist}
    while true
        print(terminal, color_string(to_qe(nml), 'b'))
        # It will continuously print until the user chooses `"no"`, i.e., he/she is satisfied.
        isdone = pairs((false, true))[request(
            terminal,
            color_string(
                "We generate an example `$(nameof(T))`. Want to change/add any field?",
                'r',
            ),
            RadioMenu(["yes", "no"]),
        )]
        if !isdone
            while true
                print(terminal, c"Type a field name: "r)
                field = strip(readline(terminal)) |> Symbol
                # Once a field successfully changes, go back to the above menu.
                # The code will asks the user whether to change another field.
                if hasfield(T, field)
                    print(terminal, c"Type its value: "r)
                    nml = set(
                        nml,
                        PropertyLens{field}(),
                        parse(fieldtype(T, field), readline(terminal)),
                    )
                    break
                end
                # If the field has a wrong name, go back to `"Type a field name: "`.
                println(terminal, c"Unknown field given! Try again!"r)
                continue
            end
            continue
        end
        break
    end
    return nml
end # function change_field_helper

module PWscf

using REPL.Terminals: TTYTerminal
using REPL.TerminalMenus: RadioMenu, request

using QuantumESPRESSO.Namelists.PWscf:
    ControlNamelist, SystemNamelist, ElectronsNamelist, IonsNamelist, CellNamelist

using ...Wizard: @c_str
using ..Namelists

function Namelists.namelist_helper(
    terminal::TTYTerminal,
    ::Type{T},
) where {T<:ControlNamelist}
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
    return Namelists.setfield_helper(terminal, control)
end # function namelist_helper
function Namelists.namelist_helper(
    terminal::TTYTerminal,
    ::Type{T},
) where {T<:SystemNamelist}
    print(terminal, c"Please input the Bravais lattice index `ibrav`: "r)
    ibrav = parse(Int, readline(terminal))
    print(terminal, c"Please input a `celldm` 1-6 (separated by spaces): "r)
    celldm = map(x -> parse(Float64, x), split(readline(terminal), " ", keepempty = false))
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
    return Namelists.setfield_helper(terminal, system)
end # function namelist_helper
function Namelists.namelist_helper(
    terminal::TTYTerminal,
    ::Type{T},
) where {T<:ElectronsNamelist}
    print(
        terminal,
        c"Please input the convergence threshold for selfconsistency `conv_thr`: "r,
    )
    conv_thr = parse(Float64, readline(terminal))
    diagonalizations = pairs(("david", "cg", "cg-serial", "david-serial"))
    diagonalization = diagonalizations[request(
        terminal,
        c"Please input the diagonalization method `diagonalization`: "r,
        RadioMenu(collect(values(diagonalizations))),
    )]
    electrons = T(conv_thr = conv_thr, diagonalization = diagonalization)
    return Namelists.setfield_helper(terminal, electrons)
end # function namelist_helper
function Namelists.namelist_helper(terminal::TTYTerminal, ::Type{T}) where {T<:IonsNamelist}
    ion_dynamics_pool =
        pairs(("none", "bfgs", "damp", "verlet", "langevin", "langevin-smc", "beeman"))
    ion_dynamics = ion_dynamics_pool[request(
        terminal,
        c"Please input the type of ionic dynamics `ion_dynamics`: "r,
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
        c"Please input the ions temperature `ion_temperature`: "r,
        RadioMenu(collect(values(ion_temperature_pool))),
    )]
    ions = T(ion_dynamics = ion_dynamics, ion_temperature = ion_temperature)
    return Namelists.setfield_helper(terminal, ions)
end # function namelist_helper
function Namelists.namelist_helper(terminal::TTYTerminal, ::Type{T}) where {T<:CellNamelist}
    cell_dynamics_pool = pairs(("none", "sd", "damp-pr", "damp-w", "bfgs", "pr", "w"))
    cell_dynamics = cell_dynamics_pool[request(
        terminal,
        c"Please input the type of dynamics for the cell `cell_dynamics`: "r,
        RadioMenu(collect(values(cell_dynamics_pool))),
    )]
    print(
        terminal,
        c"Please input the target pressure [KBar] in a variable-cell md or relaxation run `press`: "r,
    )
    press = parse(Float64, readline(terminal))
    print(
        terminal,
        c"Please input the fictitious cell mass [amu] for variable-cell simulations `wmass`: "r,
    )
    wmass = parse(Float64, readline(terminal))
    print(
        terminal,
        c"Please input the Convergence threshold on the pressure for variable cell `press_conv_thr`: "r,
    )
    press_conv_thr = parse(Float64, readline(terminal))
    cell = T(
        cell_dynamics = cell_dynamics,
        press = press,
        wmass = wmass,
        press_conv_thr = press_conv_thr,
    )
    return Namelists.setfield_helper(terminal, cell)
end # function namelist_helper

end # module PWscf

module PHonon

using REPL.Terminals: TTYTerminal
using REPL.TerminalMenus: RadioMenu, request

using QuantumESPRESSO.Namelists.PHonon:
    PhNamelist, Q2rNamelist, MatdynNamelist, DynmatNamelist

using ...Wizard: @c_str
using ..Namelists

function Namelists.namelist_helper(terminal::TTYTerminal, ::Type{T}) where {T<:PhNamelist}
    print(
        terminal,
        c"Please input the atomic mass [amu] of each atomic type `amass` (separated by spaces): "r,
    )
    amass = map(x -> parse(Float64, x), split(readline(terminal), " ", keepempty = false))
    epsil_pool = pairs((false, true))
    epsil = epsil_pool[request(
        terminal,
        c"Please select the `epsil`: "r,
        RadioMenu([false, true]),
    )]
    q_in_band_form_pool = pairs((false, true))
    q_in_band_form = q_in_band_form_pool[request(
        terminal,
        c"Please select the `q_in_band_form`: "r,
        RadioMenu([false, true]),
    )]
    print(
        terminal,
        c"Please input parameters of the Monkhorst-Pack grid `nq` 1-3 (separated by spaces): "r,
    )
    nq1, nq2, nq3 =
        map(x -> parse(Float64, x), split(readline(terminal), " ", keepempty = false))
    print(
        terminal,
        c"Please input parameters of the Monkhorst-Pack grid `nk` 1-3 (separated by spaces): "r,
    )
    nk1, nk2, nk3 =
        map(x -> parse(Float64, x), split(readline(terminal), " ", keepempty = false))
    print(terminal, c"Please input offset `k` 1-3 (separated by spaces): "r)
    k1, k2, k3 =
        map(x -> parse(Float64, x), split(readline(terminal), " ", keepempty = false))
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
    return Namelists.setfield_helper(terminal, ph)
end # function namelist_helper
function Namelists.namelist_helper(terminal::TTYTerminal, ::Type{T}) where {T<:Q2rNamelist}
    print(terminal, c"name of input dynamical matrices `fildyn`: "r)
    fildyn = strip(readline(terminal))
    print(terminal, c"name of output force constants `flfrc`: "r)
    flfrc = strip(readline(terminal))
    zasr_pool = pairs(("no", "simple", "crystal", "one-dim", "zero-dim"))
    zasr = zasr_pool[request(
        terminal,
        c"Please input the type of acoustic sum rules used for the Born effective charges `zasr`: "r,
        RadioMenu(collect(values(zasr_pool))),
    )]
    q2r = T(fildyn = fildyn, flfrc = flfrc, zasr = zasr)
    return Namelists.setfield_helper(terminal, q2r)
end # function namelist_helper
function Namelists.namelist_helper(
    terminal::TTYTerminal,
    ::Type{T},
) where {T<:MatdynNamelist}
    dos_pool = pairs((false, true))
    dos = dos_pool[request(
        terminal,
        c"Please select if calculate phonon density of states `dos`: "r,
        RadioMenu([false, true]),
    )]
    print(terminal, c"Please input the energy step, in cm^(-1) `deltaE`: "r)
    deltaE = parse(Float64, readline(terminal))
    print(
        terminal,
        c"Please input uniform q-point grid for DOS calculation `nk` 1-3 (separated by spaces): "r,
    )
    nk1, nk2, nk3 =
        map(x -> parse(Float64, x), split(readline(terminal), " ", keepempty = false))
    asr_pool = pairs(("no", "simple", "crystal", "one-dim", "zero-dim"))
    asr = asr_pool[request(
        terminal,
        c"Please input the type of acoustic sum rule `asr`: "r,
        RadioMenu(collect(values(asr_pool))),
    )]
    print(terminal, c"name of output force constants `flfrc`: "r)
    flfrc = strip(readline(terminal))
    print(terminal, c"name of input dynamical matrices `fildyn`: "r)
    fildyn = strip(readline(terminal))
    print(
        terminal,
        c"Please input the masses of atoms in the supercell (a.m.u.) `amass` (separated by spaces): "r,
    )
    amass = map(x -> parse(Float64, x), split(readline(terminal), " ", keepempty = false))
    print(terminal, c"Please input the number of atom types in the supercell `ntyp`: "r)
    ntyp = parse(Int, readline(terminal))
    q_in_band_form_pool = pairs((false, true))
    q_in_band_form = q_in_band_form_pool[request(
        terminal,
        c"Please select the `q_in_band_form`: "r,
        RadioMenu([false, true]),
    )]
    q_in_cryst_coord_pool = pairs((false, true))
    q_in_cryst_coord = q_in_cryst_coord_pool[request(
        terminal,
        c"Please select the `q_in_cryst_coord`: "r,
        RadioMenu([false, true]),
    )]
    nosym_pool = pairs((false, true))
    nosym = nosym_pool[request(
        terminal,
        c"Please select if impose symmetry and time reversal `nosym`: "r,
        RadioMenu([false, true]),
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
        nosym = nosym,
    )
    return Namelists.setfield_helper(terminal, matdyn)
end # function namelist_helper
function Namelists.namelist_helper(
    terminal::TTYTerminal,
    ::Type{T},
) where {T<:DynmatNamelist}
    asr_pool = pairs(("no", "simple", "crystal", "one-dim", "zero-dim"))
    asr = asr_pool[request(
        terminal,
        c"Please select the type of acoustic sum rule `asr`: "r,
        RadioMenu(collect(values(asr_pool))),
    )]
    print(
        terminal,
        c"Please input mass for each atom type `amass` (separated by spaces): "r,
    )
    amass = map(x -> parse(Float64, x), split(readline(terminal), " ", keepempty = false))
    dynmat = T(asr = asr, amass = amass)
    return Namelists.setfield_helper(terminal, dynmat)
end # function namelist_helper

end # module PHonon

end
