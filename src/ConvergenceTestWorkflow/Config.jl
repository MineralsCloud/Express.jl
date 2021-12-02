module Config

using AbInitioSoftwareBase.Commands: CommandConfig
using Configurations: from_dict, @option

using ...Express: Calculation, Action, UnitfulVector, myuparse

@option "ecutwfc" struct CutoffEnergy <: UnitfulVector
    values::AbstractVector
    unit::String = "Ry"
end

@option "k_mesh" struct KMesh <: UnitfulVector
    mesh::AbstractVector
    is_shift::AbstractVector = [0, 0, 0]
    function KMesh(mesh, is_shift)
        @assert all(mesh .>= 1)
        @assert all(0 <= x <= 1 for x in is_shift)
        return new(mesh, is_shift)
    end
end

@option struct Directories
    root::String = pwd()
end

@option struct FileNamePatterns
    input::String = "%s.in"
    output::String = "%s.out"
end

@option struct Save
    raw::String = "raw.json"
    status::String = ""
end

@option struct IOFiles
    dirs::Directories = Directories()
    pattern::FileNamePatterns = FileNamePatterns()
end

@option struct RuntimeConfig
    recipe::String
    template::String
    parameters::Union{CutoffEnergy,KMesh}
    files::IOFiles = IOFiles()
    save::Save = Save()
    cli::CommandConfig
    function RuntimeConfig(recipe, template, parameters, files, save, cli)
        @assert recipe in ("ecut", "k_mesh")
        if !isfile(template)
            @warn "I cannot find template file `$template`!"
        end
        return new(recipe, template, parameters, files, save, cli)
    end
end

end
