module Config

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

end
