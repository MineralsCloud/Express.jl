"""
# module BandStructure



# Examples

```jldoctest
julia>
```
"""
module BandStructure

using ConstructionBase: setproperties
using Crystallography.Symmetry: genpath
using QuantumESPRESSO.Inputs: InputFile
using QuantumESPRESSO.Inputs.PWscf: SpecialKPoint, KPointsCard, BandsNamelist, PWInput

import ..Step

function update_kpoints(template::PWInput, path::AbstractVector{<:AbstractVector})
    data = map(x -> SpecialKPoint(x, 1), path)
    return setproperties(template, k_points = KPointsCard("crystal_b", data))
end # function update_kpoints

# This is a helper function and should not be exported
which_calculation(step::Step{1}) = "nscf"
which_calculation(step::Step{2}) = "bands"

# This is a helper function and should not be exported
function set_calculation(step::Step, template::PWInput)
    @warn "We will set the calculation type to be \"$type\" for you."
    return setproperties(
        template,
        setproperties(template.control, calculation = which_calculation(step)),
    )
end # function set_calculation

function prepare(step::Step{1}, inputs::AbstractVector{<:AbstractString}, template::PWInput)
    # Checking parameters
    template = set_calculation(step, template)
end# function prepare
function prepare(
    step::Step{2},
    inputs::AbstractVector{<:AbstractString},
    template::PWInput,
    nodes::AbstractVector{<:AbstractVector},
    densities::AbstractVector{<:Integer} = 100 * ones(Int, length(nodes)),
)
    # Checking parameters
    template = set_calculation(step, template)
    template = update_kpoints(template, genpath(nodes, densities))
    for input in inputs
        write(InputFile(input), template)
    end
end# function prepare
function prepare(
    step::Step{3},
    inputs::AbstractVector{<:AbstractString},
    template::BandsNamelist,
)
    for input in inputs
        write(InputFile(input), template)
    end
end # function prepare

end
