"""
# module SelfConsistentField



# Examples

```jldoctest
julia>
```
"""
module SelfConsistentField

import JSON
using QuantumESPRESSO.Inputs.PWscf: PWscfInput

export write_metadata

function write_metadata(file::AbstractString, input::AbstractString, object::PWscfInput)
    metadata = Dict{String,Any}(
        "outdir" => object.control.outdir,
        "prefix" => object.control.prefix,
        "pseudo_dir" => object.control.pseudo_dir,
        "pseudopotentials" => [getfield(x, :pseudopotential) for x in object.atomic_species.data],
        "input" => input,
    )
    metadata["wfcdir"] = object.control.wf_collect ? metadata["outdir"] :
                         object.control.wfcdir
    metadata["wfc_namepattern"] = metadata["prefix"] * ".wfc"
    if object.control.lkpoint_dir
        metadata["lkpoint_dir"] = metadata["outdir"] * "/" * metadata["prefix"] * ".save"
    end
    lowercase(splitext(file)[2]) != ".json" && error("The file to be dumped must be a JSON file!")
    ispath(file) || touch(file)
    open(file, "r+") do io
        JSON.print(io, metadata)
    end
end # function write_metadata

end
