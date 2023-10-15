using AbInitioSoftwareBase: Input
using EasyJobsBase: Job
using ExpressBase:
    VariableCellOptimization,
    LatticeDynamics,
    SelfConsistentField,
    DensityFunctionalPerturbationTheory,
    RealSpaceForceConstants,
    PhononDispersion,
    VDOS,
    Action,
    DownloadPotentials,
    RunCmd,
    WriteInput
using ExpressBase.Files: save

struct CreateInput{T} <: Action{T}
    calculation::T
end

struct ExtractData{T} <: Action{T}
    calculation::T
end

struct GatherData{T} <: Action{T}
    calculation::T
end
(::GatherData)(iter) = collect(iter)

struct SaveData{T} <: Action{T}
    calculation::T
end
function (action::SaveData)(path, raw_data)
    raw_data = sort(collect(raw_data))  # In case the data is not sorted
    data = Dict(
        "volume" => (string ∘ first).(raw_data), "energy" => (string ∘ last).(raw_data)
    )
    return save(path, data)
end
(action::SaveData)(path) = Base.Fix1(action, path)

function parsecell end
