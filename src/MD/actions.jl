using AbInitioSoftwareBase: Input
using EasyJobsBase: Job
using ExpressBase:
    MolecularDynamics,
    IonDynamics,
    VariableCellMolecularDynamics,
    Action,
    DownloadPotentials,
    RunCmd,
    WriteInput
using ExpressBase.Files: save

struct CreateInput{T} <: Action{T}
    calculation::T
end
(action::CreateInput)(template::Input) = Base.Fix1(action, template)

struct ExtractTrajectory{T} <: Action{T}
    calculation::T
end

struct GatherTrajectories{T} <: Action{T}
    calculation::T
end
(::GatherTrajectories)(iter) = collect(iter)

struct SaveTrajectory{T} <: Action{T}
    calculation::T
end
function (action::SaveTrajectory)(path, raw)
    raw = collect(raw)  # In case the data is not sorted
    data = Dict("volume" => (string ∘ first).(raw), "energy" => (string ∘ last).(raw))
    return save(path, data)
end
(action::SaveTrajectory)(path) = Base.Fix1(action, path)
