module Recipes

using AbInitioSoftwareBase: load
using Serialization: deserialize
using SimpleWorkflows: Job, Workflow, run!, →, ⇉, ⭃
using ..PhononWorkflow:
    Scf,
    Dfpt,
    RealSpaceForceConstants,
    PhononDispersion,
    VDos,
    DownloadPotentials,
    LogMsg,
    MakeInput,
    RunCmd,
    buildjob

export buildworkflow, run!

function buildworkflow(cfgfile)
    dict = load(cfgfile)
    if isfile(dict["save"]["status"])
        w = deserialize(dict["save"]["status"])
        typeassert(w, Workflow)
        return w
    else
        x = if dict["recipe"] == "phonon dispersion"
            PhononDispersion
        elseif dict["recipe"] == "vdos"
            VDos
        else
            error("unsupported option!")
        end
        a0 = buildjob(DownloadPotentials{Scf}(), cfgfile)
        a = Job(() -> LogMsg{Scf}()(; start = true))
        b = buildjob(MakeInput{Scf}(), cfgfile)
        c = buildjob(RunCmd{Scf}(), cfgfile)
        d = Job(() -> LogMsg{Scf}()(; start = false))
        e = Job(() -> LogMsg{Dfpt}()(; start = true))
        f = buildjob(MakeInput{Dfpt}(), cfgfile)
        g = buildjob(RunCmd{Dfpt}(), cfgfile)
        h = Job(() -> LogMsg{Dfpt}()(; start = false))
        i = Job(() -> LogMsg{RealSpaceForceConstants}()(; start = true))
        j = buildjob(MakeInput{RealSpaceForceConstants}(), cfgfile)
        k = buildjob(RunCmd{RealSpaceForceConstants}(), cfgfile)
        l = Job(() -> LogMsg{RealSpaceForceConstants}()(; start = false))
        m = Job(() -> LogMsg{x}()(; start = true))
        n = buildjob(MakeInput{x}(), cfgfile)
        o = buildjob(RunCmd{x}(), cfgfile)
        p = Job(() -> LogMsg{x}()(; start = false))
        a0 → a ⇉ b → c ⭃ d → e ⇉ f → g ⭃ h → i ⇉ j → k ⭃ l → m ⇉ n → o ⭃ p
        return Workflow(a0)
    end
end

end
