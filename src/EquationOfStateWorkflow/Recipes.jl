module Recipes

using Serialization: deserialize
using SimpleWorkflows: AtomicJob, Workflow, run!, →

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
        b = AtomicJob(() -> MakeInput{Scf}()(cfgfile))
        c = AtomicJob(() -> RunCmd{Scf}()(cfgfile))
        d = AtomicJob(() -> FitEos{Scf}()(cfgfile))
        e = AtomicJob(() -> GetData{Scf}()(string(Scf) * ".json", last.(config.files)))
        f = AtomicJob(() -> LogMsg{Scf}()(false))
        g = AtomicJob(() -> LogMsg{VcOptim}()(true))
        h = AtomicJob(() -> MakeInput{VcOptim}()(cfgfile))
        i = AtomicJob(() -> RunCmd{VcOptim}()(cfgfile))
        j = AtomicJob(() -> FitEos{VcOptim}()(cfgfile))
        k = AtomicJob(
            () -> GetData{VcOptim}()(string(VcOptim) * ".json", last.(config.files)),
        )
        l = AtomicJob(() -> LogMsg{VcOptim}()(false))
        a → b → c → d → e → f → g → h → i → j → k → l
        return Workflow(a, b, c, d, e, f, g, h, i, j, k, l)
    end
end

end
