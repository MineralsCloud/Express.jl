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
using SimpleWorkflow: InternalAtomicJob, chain
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

buildjob(x::FitEos, args...) = InternalAtomicJob(() -> x(args...))
buildjob(::MakeInput{T}, cfgfile) where {T} =
    InternalAtomicJob(() -> MakeInput{T}()(cfgfile))
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

function _alert(pressures)
    if length(pressures) <= 5
        @info "pressures <= 5 may give unreliable results, consider more if possible!"
    end
    if minimum(pressures) >= zero(eltype(pressures))
        @warn "for better fitting, we need at least 1 negative pressure!"
    end
end

shortname(calc::ScfOrOptim) = shortname(typeof(calc))

function materialize_eos(config)
    selection = lowercase(config["type"])
    constructor = if selection in ("m", "murnaghan")
        Murnaghan
    elseif selection in ("bm2", "birchmurnaghan2nd", "birch-murnaghan-2")
        BirchMurnaghan2nd
    elseif selection in ("bm3", "birchmurnaghan3rd", "birch-murnaghan-3")
        BirchMurnaghan3rd
    elseif selection in ("bm4", "birchmurnaghan4th", "birch-murnaghan-4")
        BirchMurnaghan4th
    elseif selection in ("pt2", "poiriertarantola2nd", "poirier-tarantola-2")
        PoirierTarantola2nd
    elseif selection in ("pt3", "poiriertarantola3rd", "poirier-tarantola-3")
        PoirierTarantola3rd
    elseif selection in ("pt4", "poiriertarantola4th", "poirier-tarantola-4")
        PoirierTarantola4th
    elseif selection in ("v", "vinet")
        Vinet
    else
        error("unsupported eos type `\"$type\"`!")
    end
    values = map(config["parameters"]) do (v, u)
        v * uparse(u; unit_context = [Unitful, UnitfulAtomic])
    end
    return constructor(values...)
end

function checkconfig(config)
    for key in ("templates", "pressures", "workdir", "pressures")
        @assert haskey(config, key) "`\"$key\"` was not found in settings!"
    end
    if !haskey(config["pressures"], "unit")
        @info "no unit provided for `\"pressures\"`! \"GPa\" is assumed!"
    end
    @assert haskey(config, "trial_eos") || haskey(config, "volumes") "either `\"trial_eos\"` or `\"volumes\"` is required in settings!"
    if !isdir(expanduser(config["workdir"]))
        @warn "`workdir` \"$(config["workdir"])\" is not reachable, be careful!"
    end
    for path in config["templates"]
        if !isfile(path)
            @warn "template \"$path\" is not reachable, be careful!"
        end
    end
    _alert(config["pressures"]["values"])
    for key in ("type", "parameters")
        @assert haskey(config["trial_eos"], key) "the trial eos needs `\"$key\"` specified!"
    end
    checkconfig(currentsoftware(), config["qe"])  # To be implemented
end

function materialize end

module DefaultActions

using AbInitioSoftwareBase: save, load, extension
using AbInitioSoftwareBase.Inputs: Input, writeinput
using EquationsOfStateOfSolids.Collections:
    EquationOfStateOfSolids, EnergyEOS, PressureEOS, getparam
using EquationsOfStateOfSolids.Fitting: eosfit
using Serialization: serialize, deserialize

using ...Express: Action, loadconfig
using ..EosFitting: ScfOrOptim, Scf, iofiles, shortname

include("makeinput.jl")
include("getdata.jl")
include("fiteos.jl")
include("saveeos.jl")

end

using .DefaultActions: MakeInput, GetData, FitEos, SaveEos

end
