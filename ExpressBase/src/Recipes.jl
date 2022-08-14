module Recipes

import CompositionsBase: decompose, deopcompose

abstract type Recipe end
# See https://github.com/JuliaLang/julia/blob/3fa2d26/base/operators.jl#L1078-L1083 & https://github.com/JuliaGeometry/CoordinateTransformations.jl/blob/ff9ea6e/src/core.jl#L29-L32
struct ComposedRecipe{R1<:Recipe,R2<:Recipe} <: Recipe
    r1::R1
    r2::R2
end

# function build(::Type{Workflow}, r::ComposedRecipe)
#     wf1 = build(Workflow, r.r1)
#     wf2 = build(Workflow, r.r2)
#     return wf1 → wf2
# end

# See https://github.com/JuliaFunctional/CompositionsBase.jl/blob/ac505d4/src/CompositionsBase.jl#L61
decompose(r::ComposedRecipe) = (decompose(r.r1)..., decompose(r.r2)...)

# See https://github.com/JuliaFunctional/CompositionsBase.jl/blob/ac505d4/src/CompositionsBase.jl#L80
deopcompose(r::ComposedRecipe) = (deopcompose(r.r2)..., deopcompose(r.r1)...)

# See https://github.com/JuliaLang/julia/blob/3fa2d26/base/operators.jl#L1088
Base.:∘(r1::Recipe, r2::Recipe) = ComposedRecipe(r1, r2)

# See https://github.com/JuliaGeometry/CoordinateTransformations.jl/blob/ff9ea6e/src/core.jl#L34
Base.show(io::IO, r::ComposedRecipe) = print(io, '(', r.r1, " ∘ ", r.r2, ')')

end
