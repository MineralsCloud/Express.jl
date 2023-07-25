using AbInitioSoftwareBase: Input, writetxt
using ExpressBase: Action, SCF
using ExpressBase.Files: save, load
using EasyJobsBase: Job
using Unitful: ustrip

using ..Express: distribute_procs
using .Config: ExpandConfig

import ..Express: jobify

struct CreateInput{T} <: Action{T}
    calculation::T
end


function parseoutput end

struct ExtractData{T} <: Action{T}
    calculation::T
end

end

struct TestConvergence{T} <: Action{T}
    calculation::T
end
(x::TestConvergence)(data) = isconvergent(data)

function isconvergent(a::AbstractVector)
    terms = abs.(diff(a))
    x, y, z = last(terms, 3)
    return all(0 <= r < 1 for r in (y / x, z / y))
end
