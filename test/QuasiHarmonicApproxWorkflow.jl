module QuasiHarmonicApproxWorkflow

using AbInitioSoftwareBase: load
using Express.QuasiHarmonicApproxWorkflow: QuasiHarmonicApprox
using Express.QuasiHarmonicApproxWorkflow.Config: ExpandConfig
using Test

@testset "Load a configuration file: silicon" begin
    dict = load("../examples/silicon/config.toml")
    config = ExpandConfig{QuasiHarmonicApprox}()(dict)
    for (path, file) in
        zip(config, ("input", "settings.yaml", "filelist.yaml", "static", "q_points"))
        @test path == normpath(joinpath(pwd(), "..", "examples", "silicon", file))
    end
end

end
