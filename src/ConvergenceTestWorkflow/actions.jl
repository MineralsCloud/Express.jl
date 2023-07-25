using AbInitioSoftwareBase: Input
using ExpressBase: Action
using ExpressBase.Files: save, load
using EasyJobsBase: Job
using Unitful: ustrip

struct CreateInput{T} <: Action{T}
    calculation::T
end
(action::CreateInput)(template::Input) = Base.Fix1(action, template)

struct ExtractData{T} <: Action{T}
    calculation::T
end

struct SaveData{T} <: Action{T}
    calculation::T
end
(action::SaveData)(path) = Base.Fix1(action, path)

struct TestConvergence{T} <: Action{T}
    calculation::T
end
(x::TestConvergence)(data) = isconvergent(data)

function isconvergent(a::AbstractVector)
    terms = abs.(diff(a))
    x, y, z = last(terms, 3)
    return all(0 <= r < 1 for r in (y / x, z / y))
end
