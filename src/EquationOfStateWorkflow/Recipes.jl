module Recipes

using Serialization: deserialize
using SimpleWorkflows: AtomicJob, Workflow, run!, →

using ...Express: loadconfig
using ..DefaultActions: LogMsg, MakeInput, MakeCmd, GetData, FitEos
using ..EquationOfStateWorkflow: Scf, VcOptim, iofiles, shortname

function buildworkflow(cfgfile)
    config = loadconfig(cfgfile)
    if isfile(config.recover)
        w = deserialize(config.recover)
        typeassert(w, Workflow)
        return w
    else
        a = AtomicJob(() -> LogMsg{Scf}()(true))
        b = AtomicJob(() -> MakeInput{Scf}()(cfgfile))
        c = AtomicJob(() -> MakeCmd{Scf}()(cfgfile))
        d = AtomicJob(() -> FitEos{Scf}()(cfgfile))
        e = AtomicJob(
            () -> GetData{Scf}()(shortname(Scf) * ".json", last.(iofiles(Scf(), cfgfile))),
        )
        f = AtomicJob(() -> LogMsg{Scf}()(false))
        g = AtomicJob(() -> LogMsg{VcOptim}()(true))
        h = AtomicJob(() -> MakeInput{VcOptim}()(cfgfile))
        i = AtomicJob(() -> MakeCmd{VcOptim}()(cfgfile))
        j = AtomicJob(() -> FitEos{VcOptim}()(cfgfile))
        k = AtomicJob(
            () -> GetData{VcOptim}()(
                shortname(VcOptim) * ".json",
                last.(iofiles(VcOptim(), cfgfile)),
            ),
        )
        l = AtomicJob(() -> LogMsg{VcOptim}()(false))
        a → b → c → d → e → f → g → h → i → j → k → l
        return Workflow(a, b, c, d, e, f, g, h, i, j, k, l)
    end
end

end
