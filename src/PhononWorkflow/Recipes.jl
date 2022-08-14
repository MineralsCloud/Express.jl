module Recipes

using AbInitioSoftwareBase: load
using ExpressBase: Scf, Dfpt, RealSpaceForceConstants, PhononDispersion, VDos
using ExpressWorkflowMaker.Templates: DownloadPotentials, LogTime, RunCmd, jobify
using SimpleWorkflows: Job, Workflow, run!, →, ⇉, ⇶, ⭃

using ..PhononWorkflow: MakeInput

export buildworkflow, run!

function buildworkflow(cfgfile)
    dict = load(cfgfile)
    x = if dict["recipe"] == "phonon dispersion"
        PhononDispersion
    elseif dict["recipe"] == "vdos"
        VDos
    else
        error("unsupported option!")
    end
    a0 = jobify(DownloadPotentials{Scf}(), cfgfile)
    a = jobify(LogTime{Scf}())
    b = jobify(MakeInput{Scf}(), cfgfile)
    c = jobify(RunCmd{Scf}(), cfgfile)
    f = jobify(MakeInput{Dfpt}(), cfgfile)
    g = jobify(RunCmd{Dfpt}(), cfgfile)
    i = jobify(LogTime{RealSpaceForceConstants}())
    j = jobify(MakeInput{RealSpaceForceConstants}(), cfgfile)
    k = jobify(RunCmd{RealSpaceForceConstants}(), cfgfile)
    m = jobify(LogTime{x}())
    n = jobify(MakeInput{x}(), cfgfile)
    o = jobify(RunCmd{x}(), cfgfile)
    a0 → a ⇉ b ⇶ c ⇉ f ⇶ g ⭃ i ⇉ j ⇶ k ⭃ m ⇉ n ⇶ o
    return Workflow(a0)
end

end
