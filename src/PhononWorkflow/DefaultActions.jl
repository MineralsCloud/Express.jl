module DefaultActions

using AbInitioSoftwareBase.Inputs: Input, writetxt
using Dates: now, format
using Logging: with_logger, current_logger

using ...Express: Action, Calculation, LatticeDynamics, Scf
using ...EquationOfStateWorkflow: VcOptim
using ..PhononWorkflow:
    Dfpt, RealSpaceForceConstants, PhononDispersion, VDos, shortname, prevcalc, order

import ..PhononWorkflow: buildjob

@action MakeCmd

include("MakeInput.jl")
include("LogMsg.jl")

end
