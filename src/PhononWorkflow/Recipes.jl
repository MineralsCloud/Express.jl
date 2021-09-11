module Recipes

using AbInitioSoftwareBase: load
using Serialization: deserialize
using SimpleWorkflows: AtomicJob, Workflow, run!, ▷, ⋲, ⋺
using ..PhononWorkflow: Scf, Dfpt, RealSpaceForceConstants, PhononDispersion, VDos
using ..DefaultActions: LogMsg, MakeInput, RunCmd, buildjob

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
        a = AtomicJob(() -> LogMsg{Scf}()(; start = true))
        b = buildjob(MakeInput{Scf}(), cfgfile)
        c = buildjob(RunCmd{Scf}(), cfgfile)
        d = AtomicJob(() -> LogMsg{Scf}()(; start = false))
        e = AtomicJob(() -> LogMsg{Dfpt}()(; start = true))
        f = buildjob(MakeInput{Dfpt}(), cfgfile)
        g = buildjob(RunCmd{Dfpt}(), cfgfile)
        h = AtomicJob(() -> LogMsg{Dfpt}()(; start = false))
        i = AtomicJob(() -> LogMsg{RealSpaceForceConstants}()(; start = true))
        j = buildjob(MakeInput{RealSpaceForceConstants}(), cfgfile)
        k = buildjob(RunCmd{RealSpaceForceConstants}(), cfgfile)
        l = AtomicJob(() -> LogMsg{RealSpaceForceConstants}()(; start = false))
        m = AtomicJob(() -> LogMsg{x}()(; start = true))
        n = buildjob(MakeInput{x}(), cfgfile)
        o = buildjob(RunCmd{x}(), cfgfile)
        p = AtomicJob(() -> LogMsg{x}()(; start = false))
        (
            (
                ((((((((((((a ⋲ b) ▷ c) ⋺ d) ▷ e) ⋲ f) ▷ g) ⋺ h) ▷ i) ⋲ j) ▷ k) ⋺ l) ▷ m) ⋲
                n
            ) ▷ o
        ) ⋺ p
        return Workflow(
            a,
            b...,
            c...,
            d,
            e,
            f...,
            g...,
            h,
            i,
            j...,
            k...,
            l,
            m,
            n...,
            o...,
            p,
        )
    end
end

end
