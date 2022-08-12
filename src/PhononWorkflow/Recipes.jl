module Recipes

using AbInitioSoftwareBase: load
using ExpressBase: Scf, Dfpt, RealSpaceForceConstants, PhononDispersion, VDos
using ExpressWorkflowMaker.Templates: DownloadPotentials, LogTime, RunCmd, jobify
using SimpleWorkflows: Job, Workflow, run!, →, ⇉, ⇶, ⭃

using ..PhononWorkflow: MakeInput, buildjob

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
    a0 = buildjob(DownloadPotentials{Scf}(), cfgfile)
    a = Job(() -> LogTime{Scf}()())
    b = buildjob(MakeInput{Scf}(), cfgfile)
    c = buildjob(RunCmd{Scf}(), cfgfile)
    d = Job(() -> LogTime{Scf}()())
    f = buildjob(MakeInput{Dfpt}(), cfgfile)
    g = buildjob(RunCmd{Dfpt}(), cfgfile)
    i = Job(() -> LogTime{RealSpaceForceConstants}()())
    j = buildjob(MakeInput{RealSpaceForceConstants}(), cfgfile)
    k = buildjob(RunCmd{RealSpaceForceConstants}(), cfgfile)
    m = Job(() -> LogTime{x}()())
    n = buildjob(MakeInput{x}(), cfgfile)
    o = buildjob(RunCmd{x}(), cfgfile)
    a0 → a ⇉ b ⇶ c ⭃ d ⇉ f ⇶ g ⭃ i ⇉ j ⇶ k ⭃ m ⇉ n ⇶ o
    return Workflow(a0)
end

end
