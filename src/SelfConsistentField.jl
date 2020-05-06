"""
# module SelfConsistentField



# Examples

```jldoctest
julia>
```
"""
module SelfConsistentField

using Parameters: @with_kw
using QuantumESPRESSO.Inputs: InputFile
using QuantumESPRESSO.Inputs.PWscf: PWInput, KPointsCard
using Setfield: set, @set!

using Express:
    ScfCalculation,
    PhononCalculation,
    StructureOptimization,
    CPMD,
    PrepareInput,
    LaunchJob,
    AnalyseOutput

Step(::ScfCalculation, ::PrepareInput) = Step(1)
Step(::ScfCalculation, ::LaunchJob) = Step(2)
Step(::ScfCalculation, ::AnalyseOutput) = Step(3)

@with_kw struct Settings
    template::String
    ecutwfc
    ecutrho
    mesh
end

function (step::Step{1})(inputs, template::PWInput, settings::Settings)
    @set! template.control.calculation = "scf"
    @set! template.control.verbosity = "high"
    @set! template.control.tstress = true
    @set! template.control.tprnfor = true
    objects = similar(inputs, PWInput)  # Create an array of `undef` of `PWInput` type
    for (i, (input, wfc, rho, mesh)) in enumerate(zip(inputs, settings.ecutwfc, settings.ecutrho, settings.mesh))
        object = @set template.system.ecutwfc = wfc
        @set! object.system.ecutrho = rho
        @set! object.k_points = KPointsCard(mesh)
        objects[i] = object  # `write` will create a file if it doesn't exist.
        mkpath(joinpath(splitpath(input)[1:end-1]...))
        touch(input)
        write(InputFile(input), object)  # Write the `object` to a Quantum ESPRESSO input file
    end
    return objects
end

end
