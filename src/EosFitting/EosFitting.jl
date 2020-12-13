module EosFitting

using Crystallography: cellvolume
using EquationsOfStateOfSolids.Collections:
    Murnaghan,
    BirchMurnaghan2nd,
    BirchMurnaghan3rd,
    BirchMurnaghan4th,
    PoirierTarantola2nd,
    PoirierTarantola3rd,
    PoirierTarantola4th,
    Vinet
using SimpleWorkflow: chain
using Unitful: uparse, ustrip, @u_str
import Unitful
import UnitfulAtomic

import ..Express
using ..Express:
    Optimization,
    SelfConsistentField,
    Scf,
    FixedIonSelfConsistentField,
    MakeCmd,
    calculation,
    currentsoftware,
    makescript,
    loadconfig

export SelfConsistentField,
    Scf,
    FixedIonSelfConsistentField,
    StructuralOptimization,
    FixedCellOptimization,
    VariableCellOptimization,
    StOptim,
    VcOptim,
    MakeInput,
    FitEos,
    GetData,
    SaveEos,
    loadconfig,
    iofiles,
    calculation,
    makescript,
    buildjob

println("Load the `EosFitting` module of the corresponding software before running some functions!")

const UNIT_CONTEXT = [Unitful, UnitfulAtomic]

struct StructuralOptimization <: Optimization end
struct VariableCellOptimization <: Optimization end
# See https://www.quantum-espresso.org/Doc/pw_user_guide/node10.html
const FixedCellOptimization = StructuralOptimization
const StOptim = StructuralOptimization
const VcOptim = VariableCellOptimization
const ScfOrOptim = Union{SelfConsistentField,Optimization}

function iofiles(T::ScfOrOptim, cfgfile)
    config = loadconfig(cfgfile)
    return map(config.dirs) do dir
        prefix = joinpath(dir, shortname(T))
        prefix * ".in" => prefix * ".out"
    end
end

function buildjob(x::MakeCmd{T}, cfgfile) where {T}
    config = loadconfig(cfgfile)
    io = iofiles(T(), cfgfile)
    infiles, outfiles = first.(io), last.(io)
    return Express.buildjob(
        x,
        outfiles,
        infiles,
        config.manager.np,
        config.bin;
        use_shell = config.use_shell,
    )
end

function buildworkflow(cfgfile)
    step1 = buildjob(MakeInput{Scf}(), cfgfile)
    step12 = chain(step1, buildjob(MakeCmd{Scf}(), cfgfile)[1])
    step123 = chain(step12[end], buildjob(FitEos{Scf}(), cfgfile))
    step4 = buildjob(MakeInput{VcOptim}(), cfgfile)
    step45 = chain(step4, buildjob(MakeCmd{VcOptim}(), cfgfile)[1])
    step456 = chain(step45[end], buildjob(FitEos{VcOptim}(), cfgfile))
    step16 = chain(step123[end], step456[1])
    return step16
end

shortname(calc::ScfOrOptim) = shortname(typeof(calc))

function materialize_eos(config)
    name = lowercase(config["name"])
    ctor = if name in ("m", "murnaghan")
        Murnaghan
    elseif name in ("bm2", "birchmurnaghan2nd", "birch-murnaghan-2")
        BirchMurnaghan2nd
    elseif name in ("bm3", "birchmurnaghan3rd", "birch-murnaghan-3")
        BirchMurnaghan3rd
    elseif name in ("bm4", "birchmurnaghan4th", "birch-murnaghan-4")
        BirchMurnaghan4th
    elseif name in ("pt2", "poiriertarantola2nd", "poirier-tarantola-2")
        PoirierTarantola2nd
    elseif name in ("pt3", "poiriertarantola3rd", "poirier-tarantola-3")
        PoirierTarantola3rd
    elseif name in ("pt4", "poiriertarantola4th", "poirier-tarantola-4")
        PoirierTarantola4th
    elseif name in ("v", "vinet")
        Vinet
    else
        error("unsupported eos type `\"$type\"`!")
    end
    values = (
        v * uparse(string(u); unit_context = UNIT_CONTEXT) for
        (v, u) in config["parameters"]
    )
    return ctor(values...)
end

function materialize_press(config)
    unit = uparse(
        if haskey(config, "unit")
            config["unit"]
        else
            @info "no unit provided for `\"pressures\"`! \"GPa\" is assumed!"
            u"GPa"
        end;
        unit_context = UNIT_CONTEXT,
    )
    return map(Base.Fix1(*, unit), config["values"])
end

function materialize_vol(config, templates)
    if haskey(config, "volumes")
        subconfig = config["volumes"]
        unit = uparse(
            if haskey(subconfig, "unit")
                subconfig["unit"]
            else
                @info "no unit provided for `\"volumes\"`! \"bohr^3\" is assumed!"
                u"bohr^3"
            end;
            unit_context = UNIT_CONTEXT,
        )
        return map(Base.Fix1(*, unit), subconfig["values"])
    else
        return map(cellvolume, templates) * u"bohr^3"
    end
end

function materialize_dirs(config, pressures)
    return map(pressures) do pressure
        abspath(joinpath(expanduser(config), "p" * string(ustrip(pressure))))
    end
end

function _alert(pressures)
    if length(pressures) <= 5
        @info "pressures <= 5 may give unreliable results, consider more if possible!"
    end
    if minimum(pressures) >= zero(eltype(pressures))
        @warn "for better fitting, we need at least 1 negative pressure!"
    end
end

function checkconfig(config)
    for key in ("pressures", "qe", "templates", "workdir")
        @assert haskey(config, key) "`\"$key\"` was not found in config!"
    end
    checkconfig(currentsoftware(), config["qe"])  # To be implemented
    let subconfig = config["pressures"], values = subconfig["values"]
        _alert(values)
        if length(config["templates"]) != 1
            if length(values) != length(config["templates"])
                throw(DimensionMismatch("templates and pressures have different lengths!"))
            end
        end
    end
    let workdir = expanduser(config["workdir"])
        if !isdir(workdir)
            @warn "`\"workdir\"` \"$workdir\" is not reachable, be careful!"
        end
    end
    for path in config["templates"]
        if !isfile(path)
            @warn "template \"$path\" is not reachable, be careful!"
        end
    end
    if haskey(config, "trial_eos")
        @assert !haskey(config, "volumes") "key \"trial_eos\" and \"volumes\" are mutually exclusive!"
        for key in ("name", "parameters")
            @assert haskey(config["trial_eos"], key) "the trial eos needs `\"$key\"` specified!"
        end
    end
    if haskey(config, "volumes")
        subconfig = config["volumes"]
        if length(subconfig["values"]) != 1
            if length(subconfig["values"]) != length(config["pressures"]["values"])
                throw(DimensionMismatch("volumes and pressures have different lengths!"))
            end
        end
    end
    return
end

function materialize end

module DefaultActions

using AbInitioSoftwareBase: save, load, extension
using AbInitioSoftwareBase.Inputs: Input, writeinput
using EquationsOfStateOfSolids.Collections:
    EquationOfStateOfSolids, EnergyEOS, PressureEOS, Parameters, getparam
using EquationsOfStateOfSolids.Fitting: eosfit
using Serialization: serialize, deserialize
using SimpleWorkflow: InternalAtomicJob

using ...Express: Action, loadconfig
using ..EosFitting: ScfOrOptim, Scf, iofiles, shortname
import ...EosFitting: buildjob

include("makeinput.jl")
include("getdata.jl")
include("fiteos.jl")
include("saveeos.jl")

end

using .DefaultActions: MakeInput, GetData, FitEos, SaveEos

end
