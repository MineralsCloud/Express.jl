module EosFitting

using AbInitioSoftwareBase: save, load, extension
using AbInitioSoftwareBase.CLI: Mpiexec
using AbInitioSoftwareBase.Inputs: Input, writeinput
using Compat: isnothing
using EquationsOfStateOfSolids.Collections:
    EquationOfStateOfSolids,
    EnergyEOS,
    PressureEOS,
    Parameters,
    Murnaghan,
    BirchMurnaghan2nd,
    BirchMurnaghan3rd,
    BirchMurnaghan4th,
    PoirierTarantola2nd,
    PoirierTarantola3rd,
    PoirierTarantola4th,
    Vinet,
    getparam
using EquationsOfStateOfSolids.Fitting: eosfit
using Serialization: serialize, deserialize
using SimpleWorkflow: InternalAtomicJob, chain
using Unitful: uparse
import Unitful
import UnitfulAtomic

import ..Express
using ..Express:
    Calculation,
    Optimization,
    SelfConsistentField,
    Scf,
    FixedIonSelfConsistentField,
    Action,
    MakeCmd,
    calculation,
    distprocs,
    makescript,
    load_settings

export SelfConsistentField,
    Scf,
    FixedIonSelfConsistentField,
    StructuralOptimization,
    FixedCellOptimization,
    VariableCellOptimization,
    StOptim,
    VcOptim,
    load_settings,
    makeinput,
    fiteos,
    saveeos,
    iofiles,
    getdata,
    calculation,
    makescript,
    writeinput,
    buildjob

struct StructuralOptimization <: Optimization end
struct VariableCellOptimization <: Optimization end
# See https://www.quantum-espresso.org/Doc/pw_user_guide/node10.html
const FixedCellOptimization = StructuralOptimization
const StOptim = StructuralOptimization
const VcOptim = VariableCellOptimization
const ScfOrOptim = Union{SelfConsistentField,Optimization}

include("makeinput.jl")
include("getdata.jl")
include("fiteos.jl")
include("saveeos.jl")

function iofiles(T::ScfOrOptim, cfgfile)
    settings = load_settings(cfgfile)
    return map(settings.dirs) do dir
        prefix = joinpath(dir, shortname(T))
        prefix * ".in" => prefix * ".out"
    end
end

buildjob(x::FitEos, args...) = InternalAtomicJob(() -> x(args...))
buildjob(::MakeInput{T}, cfgfile) where {T} =
    InternalAtomicJob(() -> MakeInput(T())(cfgfile))
function buildjob(x::MakeCmd{T}, cfgfile) where {T}
    settings = load_settings(cfgfile)
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
    step1 = buildjob(MakeInput(SelfConsistentField()), cfgfile)
    step12 = chain(step1, buildjob(MakeCmd(SelfConsistentField()), cfgfile)[1])
    step123 = chain(step12[end], buildjob(FitEos(SelfConsistentField()), cfgfile))
    step4 = buildjob(MakeInput(VariableCellOptimization()), cfgfile)
    step45 = chain(step4, buildjob(MakeCmd(VariableCellOptimization()), cfgfile)[1])
    step456 = chain(step45[end], buildjob(FitEos(VariableCellOptimization()), cfgfile))
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

function expand_settings end

function check_software_settings end

shortname(calc::ScfOrOptim) = shortname(typeof(calc))

vscaling()::NTuple{2,<:AbstractFloat} = (0.5, 1.5)

function expandeos(settings)
    type = lowercase(settings["type"])
    constructor = if type in ("m", "murnaghan")
        Murnaghan
    elseif type in ("bm2", "birchmurnaghan2nd", "birch-murnaghan-2")
        BirchMurnaghan2nd
    elseif type in ("bm3", "birchmurnaghan3rd", "birch-murnaghan-3")
        BirchMurnaghan3rd
    elseif type in ("bm4", "birchmurnaghan4th", "birch-murnaghan-4")
        BirchMurnaghan4th
    elseif type in ("pt2", "poiriertarantola2nd", "poirier-tarantola-2")
        PoirierTarantola2nd
    elseif type in ("pt3", "poiriertarantola3rd", "poirier-tarantola-3")
        PoirierTarantola3rd
    elseif type in ("pt4", "poiriertarantola4th", "poirier-tarantola-4")
        PoirierTarantola4th
    elseif type in ("v", "vinet")
        Vinet
    end
    values = map(settings["parameters"]) do (v, u)
        v * uparse(u; unit_context = [Unitful, UnitfulAtomic])
    end
    return constructor(values...)
end

function check_settings(settings)
    for key in ("templates", "pressures", "trial_eos", "workdir")
        @assert haskey(settings, key)
    end
    @assert haskey(settings, "pressures") || haskey(settings, "volumes")
    if !isdir(expanduser(settings["workdir"]))
        @warn "`workdir` is not reachable, be careful!"
    end
    for path in settings["templates"]
        if !isfile(path)
            @warn "template \"$path\" is not reachable, be careful!"
        end
    end
    _alert(settings["pressures"]["values"])
    for key in ("type", "parameters")
        @assert haskey(settings["trial_eos"], key)
    end
end

end
