module Express

export Step

struct Step{N} end
Step(n) = Step{n}()

include("StructureOptimization.jl")

end # module
