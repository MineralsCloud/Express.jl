module QuasiHarmonicApprox

using ExpressBase.Files: load
using Express.QuasiHarmonicApprox: QuasiHarmonicApproximation
using Express.QuasiHarmonicApprox.Config: ExpandConfig
using Test

@testset "Load a configuration file: silicon" begin
    dict = load("../examples/silicon/config.toml")
    config = ExpandConfig{QuasiHarmonicApproximation}()(dict)
    for (path, file) in
        zip(config, ("input", "settings.yaml", "filelist.yaml", "static", "q_points"))
        @test path == normpath(joinpath(pwd(), "..", "examples", "silicon", file))
    end
end

end
