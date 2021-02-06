module EosFitting

using SimpleWorkflow: ExternalAtomicJob, chain, parallel, run!

import ..Express
using ..Express:
    Optimization,
    SelfConsistentField,
    Scf,
    FixedIonSelfConsistentField,
    calculation,
    currentsoftware,
    makescript,
    loadconfig
using ..Shell: MakeCmd

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
    run!,
    buildjob

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

function buildjob end

function buildworkflow(cfgfile)
    step1 = buildjob(MakeInput{Scf}(), cfgfile)
    step12 = chain(step1, buildjob(MakeCmd{Scf}(), cfgfile)[1])
    step123 = chain(step12[end], buildjob(FitEos{Scf}(), cfgfile))
    step3p =
        buildjob(GetData{Scf}(), shortname(Scf) * ".json", last.(iofiles(Scf(), cfgfile)))
    step123 = chain(step123[end], step3p)
    step4 = buildjob(MakeInput{VcOptim}(), cfgfile)
    step45 = chain(step4, buildjob(MakeCmd{VcOptim}(), cfgfile)[1])
    step456 = chain(step45[end], buildjob(FitEos{VcOptim}(), cfgfile))
    step6p = buildjob(
        GetData{VcOptim}(),
        shortname(VcOptim) * ".json",
        last.(iofiles(VcOptim(), cfgfile)),
    )
    step456 = chain(step456[end], step6p)
    step16 = chain(step123[end], step456[1])
    return step16
end

shortname(calc::ScfOrOptim) = shortname(typeof(calc))

include("Config.jl")

module DefaultActions

using AbInitioSoftwareBase: save, load, extension
using AbInitioSoftwareBase.Inputs: Input, writetxt
using EquationsOfStateOfSolids:
    EquationOfStateOfSolids, EnergyEquation, PressureEquation, Parameters, getparam
using EquationsOfStateOfSolids.Fitting: eosfit
using Serialization: serialize, deserialize
using SimpleWorkflow: InternalAtomicJob
using Unitful: ustrip, unit

using ...Express: Action, loadconfig
using ..EosFitting: ScfOrOptim, Scf, iofiles, shortname
import ...EosFitting: buildjob

include("MakeInput.jl")
include("GetData.jl")
include("FitEos.jl")
include("SaveEos.jl")

end

using .DefaultActions: MakeInput, GetData, FitEos, SaveEos

end
