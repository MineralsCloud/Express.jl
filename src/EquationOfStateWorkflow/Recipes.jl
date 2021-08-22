module Recipes

using AbInitioSoftwareBase: load
using Serialization: deserialize
using SimpleWorkflows: AtomicJob, Workflow, run!, ▷, ⋲, ⋺

using ..Config: ExpandConfig
using ..DefaultActions: LogMsg, MakeInput, RunCmd, GetData, FitEos, buildjob
using ..EquationOfStateWorkflow: Scf, VcOptim

export buildworkflow

function buildworkflow(cfgfile)
    dict = load(cfgfile)
    if isfile(dict["save"]["status"])
        w = deserialize(dict["save"]["status"])
        typeassert(w, Workflow)
        return w
    else
        a = AtomicJob(() -> LogMsg{Scf}()(true))
        b = buildjob(MakeInput{Scf}(), cfgfile)
        c = buildjob(RunCmd{Scf}(), cfgfile)
        d = buildjob(FitEos{Scf}(), cfgfile)
        f = AtomicJob(() -> LogMsg{Scf}()(false))
        g = AtomicJob(() -> LogMsg{VcOptim}()(true))
        h = buildjob(MakeInput{VcOptim}(), cfgfile)
        i = buildjob(RunCmd{VcOptim}(), cfgfile)
        j = buildjob(FitEos{Scf}(), cfgfile)
        l = AtomicJob(() -> LogMsg{VcOptim}()(false))
        ⋲(a, b...)
        for (x, y) in zip(b, c)
            x ▷ y
        end
        ⋺(c..., d)
        d ▷ f ▷ g
        ⋲(g, h...)
        for (x, y) in zip(h, i)
            x ▷ y
        end
        ⋺(i..., j)
        j ▷ l
        return Workflow(a, b..., c..., d, f, g, h..., i..., j, l)
    end
end

end
