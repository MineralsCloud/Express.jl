module Recipes

using CompositionsBase: ⨟
using ExpressBase.Recipes: Recipe
using Test: @testset, @test

struct RecipeA <: Recipe end
struct RecipeB <: Recipe end
struct RecipeC <: Recipe end
struct RecipeD <: Recipe end

const a = RecipeA()
const b = RecipeB()
const c = RecipeC()
const d = RecipeD()

@testset "Test `==` of `ComposedRecipe`" begin
    @test a ∘ b == a ∘ b == ∘(a, b)
    @test c ∘ d != d ∘ c != c != d
    @test a ∘ c != a ∘ b != b ∘ c
    @test a ∘ b ∘ c == ∘(a, b, c) == (a ∘ b) ∘ c != a ∘ (b ∘ c)
    @test a ∘ b ∘ c ∘ d == ∘(a, b, c, d) == ((a ∘ b) ∘ c) ∘ d
    @test (a ∘ b) ∘ (c ∘ d) != ((a ∘ b) ∘ c) ∘ d != a ∘ (b ∘ c) ∘ d
end

@testset "Test `⨟` of `ComposedRecipe`" begin
    @test ⨟(a) == a
    @test a ⨟ b == ⨟(a, b) == b ∘ a
    @test ⨟(a, b, c) == a ⨟ (b ⨟ c) == (c ∘ b) ∘ a == c ∘ b ∘ a
    @test a ⨟ b ⨟ c == c ∘ (b ∘ a) != c ∘ b ∘ a
    @test ⨟(a, b, c, d) == a ⨟ (b ⨟ (c ⨟ d)) == ((d ∘ c) ∘ b) ∘ a == d ∘ c ∘ b ∘ a
    @test a ⨟ b ⨟ c ⨟ d == d ∘ (c ∘ (b ∘ a)) != d ∘ c ∘ b ∘ a
    @test (a ⨟ b) ⨟ (c ⨟ d) == (d ∘ c) ∘ (b ∘ a)
end

end
