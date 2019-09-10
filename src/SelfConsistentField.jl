"""
# module SelfConsistentField



# Examples

```jldoctest
julia>
```
"""
module SelfConsistentField

export write_metadata

function write_metadata(output::AbstractString, object::PWscfInput, input::AbstractString)
    metadata = Dict(
        "outdir" => object.control.outdir,
        "prefix" => object.control.prefix,
        "pseudo_dir" => object.control.pseudo_dir,
        "pseudopotentials" => [getfield(x, :pseudopotential) for x in object.atomic_species.data],
        "input" => input
    )
    metadata["wfcdir"] = object.control.wf_collect ? metadata["outdir"] : object.control.wfcdir
    metadata["wfc_namepattern"] = metadata["prefix"] * ".wfc"
    if object.control.lkpoint_dir
        metadata["lkpoint_dir"] = metadata["outdir"] * "/" * metadata["prefix"] * ".save"
    end
    lowercase(splitext(output)[2]) != ".json" && error("The file to be dumped must be a JSON file!")
    open(output, "r+") do io
        JSON.print(io, metadata)
    end
end # function write_metadata

end
