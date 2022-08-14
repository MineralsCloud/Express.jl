module Recipes

using CompositionsBase: ⨟, decompose, deopcompose
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

@testset "Test `decompose` of `ComposedRecipe`" begin
    @test decompose(a) == (a,)
    @test decompose(a ∘ b) == (a, b)
    @test decompose(a ∘ b ∘ c) == (a, b, c)
    @test decompose(a ∘ (b ∘ c)) == (a, b, c)
    @test decompose(a ∘ b ∘ c ∘ d) == (a, b, c, d)
    @test decompose(a ∘ (b ∘ (c ∘ d))) == (a, b, c, d)
    @test decompose((a ∘ b) ∘ (c ∘ d)) == (a, b, c, d)
    @test decompose(a ⨟ b) == decompose(⨟(a, b)) == (b, a)
    @test decompose(a ⨟ b ⨟ c) == (c, b, a)
    @test decompose(⨟(a, b, c)) == decompose(a ⨟ (b ⨟ c)) == (c, b, a)
    @test decompose(a ⨟ b ⨟ c ⨟ d) == (d, c, b, a)
    @test decompose(⨟(a, b, c, d)) == decompose(a ⨟ (b ⨟ (c ⨟ d))) == (d, c, b, a)
    @test decompose((a ⨟ b) ⨟ (c ⨟ d)) == (d, c, b, a)
end

@testset "Test `deopcompose` of `ComposedRecipe`" begin
    @test deopcompose(a) == (a,)
    @test deopcompose(a ∘ b) == (b, a)
    @test deopcompose(a ∘ b ∘ c) == (c, b, a)
    @test deopcompose(a ∘ (b ∘ c)) == (c, b, a)
    @test deopcompose(a ∘ b ∘ c ∘ d) == (d, c, b, a)
    @test deopcompose(a ∘ (b ∘ (c ∘ d))) == (d, c, b, a)
    @test deopcompose((a ∘ b) ∘ (c ∘ d)) == (d, c, b, a)
    @test deopcompose(a ⨟ b) == deopcompose(⨟(a, b)) == (a, b)
    @test deopcompose(a ⨟ b ⨟ c) == (a, b, c)
    @test deopcompose(⨟(a, b, c)) == deopcompose(a ⨟ (b ⨟ c)) == (a, b, c)
    @test deopcompose(a ⨟ b ⨟ c ⨟ d) == (a, b, c, d)
    @test deopcompose(⨟(a, b, c, d)) == deopcompose(a ⨟ (b ⨟ (c ⨟ d))) == (a, b, c, d)
    @test deopcompose((a ⨟ b) ⨟ (c ⨟ d)) == (a, b, c, d)
end

end
