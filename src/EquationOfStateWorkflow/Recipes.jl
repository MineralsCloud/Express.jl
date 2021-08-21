module Recipes

using Serialization: deserialize
using SimpleWorkflows: AtomicJob, Workflow, run!, ▷, ⋲, ⋺

using ...Config: loadconfig
using ..DefaultActions: LogMsg, MakeInput, RunCmd, GetData, FitEos
using ..EquationOfStateWorkflow: Scf, VcOptim

export buildworkflow

function buildworkflow(cfgfile)
    config = loadconfig(cfgfile)
    if isfile(config.recover)
        w = deserialize(config.recover)
        typeassert(w, Workflow)
        return w
    else
        a = AtomicJob(() -> LogMsg{Scf}()(true))
        b = map(x -> AtomicJob(() -> x), MakeInput{Scf}()(cfgfile))
        c = map(x -> AtomicJob(() -> x), RunCmd{Scf}()(cfgfile))
        d = AtomicJob(() -> FitEos{Scf}()(cfgfile))
        e = AtomicJob(() -> GetData{Scf}()(string(Scf) * ".json", last.(config.files)))
        f = AtomicJob(() -> LogMsg{Scf}()(false))
        g = AtomicJob(() -> LogMsg{VcOptim}()(true))
        h = map(x -> AtomicJob(() -> x), MakeInput{VcOptim}()(cfgfile))
        i = map(x -> AtomicJob(() -> x), RunCmd{VcOptim}()(cfgfile))
        j = AtomicJob(() -> FitEos{VcOptim}()(cfgfile))
        k = AtomicJob(
            () -> GetData{VcOptim}()(string(VcOptim) * ".json", last.(config.files)),
        )
        l = AtomicJob(() -> LogMsg{VcOptim}()(false))
        ⋲(a, b...)
        for (x, y) in zip(b, c)
            x ▷ y
        end
        ⋺(c..., d)
        d ▷ e ▷ f ▷ g
        ⋲(g, h...)
        for (x, y) in zip(h, i)
            x ▷ y
        end
        ⋺(i..., j)
        j ▷ k ▷ l
        return Workflow(a, b, c, d, e, f, g, h, i, j, k, l)
    end
end

end
