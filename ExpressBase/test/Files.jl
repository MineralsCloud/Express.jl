module Files

using ExpressBase.Files: UnsupportedExtensionError, save, load, of_format, extension
using Test: @testset, @test, @test_throws

@testset "`load` and `save`" begin
    dict = Dict(
        "a" => 1,
        "b" => 1.0e4,
        "c" => [1, 2, "3", "4.0"],
        "4" => true,
        "5.0" => "Hello",
        "dict" => Dict("4" => 1, "as" => 2),
    )
    @testset "Save to YAML" begin
        file = "test.yaml"
        save(file, dict)
        @test load(file) == dict
        file = "test.yml"
        save(file, dict)
        @test load(file) == dict
        save("mixed.YaMl", dict)
        @test load("mixed.YaMl") == dict
    end
    @testset "Save to JSON" begin
        file = "test.json"
        save(file, dict)
        @test load(file) == dict
        save("mixed.Json", dict)
        @test load("mixed.Json") == dict
    end
    @testset "Save to TOML" begin
        file = "test.toml"
        save(file, dict)
        @test load(file) == dict
        save("mixed.ToMl", dict)
        @test load("mixed.ToMl") == dict
    end
    @testset "Other extensions" begin
        @test_throws UnsupportedExtensionError save("x.aml", dict)
        @test_throws UnsupportedExtensionError save("x.jon", dict)
        @test_throws UnsupportedExtensionError save("x.Tml", dict)
    end
    @testset "Format conversions" begin
        of_format("testconvert.json", "test.toml")
        @test load("testconvert.json") == load("test.toml")
        of_format("testconvert.toml", "test.json")
        @test load("testconvert.toml") == load("test.json")
    end
end

@testset "Test files with multiple extensions" begin
    @test extension("a.b_c.d.yaml") == "yaml"
    @test extension("a.b_c.d.YaMl") == "yaml"
end

end
