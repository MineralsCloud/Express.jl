module Config

using Configurations: from_dict, @option

using ...Express: Calculation, Action, UnitfulVector, myuparse

@option "ecutwfc" struct CutoffEnergy <: UnitfulVector
    values::AbstractVector
    unit::String = "Ry"
end



end
