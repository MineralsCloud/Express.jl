module EosFitting

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
using Unitful: uparse
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
    settings = loadconfig(cfgfile)
    return map(settings.dirs) do dir
        prefix = joinpath(dir, shortname(T))
        prefix * ".in" => prefix * ".out"
    end
end

function buildjob(x::MakeCmd{T}, cfgfile) where {T}
    settings = loadconfig(cfgfile)
    io = iofiles(T(), cfgfile)
    infiles, outfiles = first.(io), last.(io)
    return Express.buildjob(
        x,
        outfiles,
        infiles,
        settings.manager.np,
        settings.bin;
        use_shell = settings.use_shell,
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
    values = (v * uparse(u; unit_context = UNIT_CONTEXT) for (v, u) in config["parameters"])
    return ctor(values...)
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
        @assert haskey(config, key) "`\"$key\"` was not found in settings!"
    end
    checkconfig(currentsoftware(), config["qe"])  # To be implemented
    let subconfig = config["pressures"]
        if !haskey(subconfig, "unit")
            @info "no unit provided for `\"pressures\"`! \"GPa\" is assumed!"
        end
        _alert(subconfig["values"])
        if config["templates"] isa Vector
            if length(subconfig["values"]) != length(config["templates"])
                throw(DimensionMismatch("templates and pressures have different length!"))
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
        for key in ("type", "parameters")
            @assert haskey(config["trial_eos"], key) "the trial eos needs `\"$key\"` specified!"
        end
    end
    if haskey(config, "volumes")
        subconfig = config["volumes"]
        if !haskey(subconfig, "unit")
            @info "no unit provided for `\"volumes\"`! \"bohr^3\" is assumed!"
        end
        if subconfig["values"] isa Vector
            if length(subconfig["values"]) != length(config["pressures"]["values"])
                throw(DimensionMismatch("volumes and pressures have different length!"))
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
