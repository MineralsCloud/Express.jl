using AbInitioSoftwareBase: Input
using ExpressBase:
    Calculation, Action, SelfConsistentField, DownloadPotentials, RunCmd, WriteInput
using ExpressBase.Files: save, load
using EasyJobsBase: Job
using Unitful: ustrip

struct CreateInput{T} <: Action{T}
    calculation::T
end

struct ExtractData{T} <: Action{T}
    calculation::T
end

struct SaveData{T} <: Action{T}
    calculation::T
end
function (action::SaveData)(path, raw_data)
    raw_data = sort(collect(raw_data))  # In case the data is not sorted
    data = Dict(
        "ecut_kmesh" => (string ∘ first).(raw_data), "energy" => (string ∘ last).(raw_data)
    )
    return save(path, data)
end
(action::SaveData)(path) = Base.Fix1(action, path)

struct TestConvergence{T} <: Action{T}
    calculation::T
end
(::TestConvergence)(threshold) = Base.Fix2(isconvergent, threshold)

function isconvergent(iter, threshold)
    last3 = last(sort(iter; by=first), 3)  # Sort a `Set` of `Pair`s by the keys, i.e., increasing cutoff energies
    min, max = extrema(last.(last3))  # Get the minimum and maximum energies from the last 3 pairs
    range = abs(max - min)
    return zero(range) <= range <= threshold
end
