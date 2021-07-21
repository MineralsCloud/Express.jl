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
        return begin
            AtomicJob(() -> LogMsg{Scf}()(true)) →
            AtomicJob(() -> MakeInput{Scf}()(cfgfile)) →
            AtomicJob(() -> MakeCmd{Scf}(), cfgfile) →
            AtomicJob(() -> FitEos{Scf}()(cfgfile)) →
            AtomicJob(
                () -> GetData{Scf}()(
                    shortname(Scf) * ".json",
                    last.(iofiles(Scf(), cfgfile)),
                ),
            ) →
            AtomicJob(() -> LogMsg{Scf}()(false)) →
            AtomicJob(() -> LogMsg{VcOptim}()(true)) →
            AtomicJob(() -> MakeInput{VcOptim}()(cfgfile)) →
            AtomicJob(() -> MakeCmd{VcOptim}(), cfgfile) →
            AtomicJob(() -> FitEos{VcOptim}()(cfgfile)) →
            AtomicJob(
                () -> GetData{VcOptim}()(
                    shortname(VcOptim) * ".json",
                    last.(iofiles(VcOptim(), cfgfile)),
                ),
            ) → AtomicJob(() -> LogMsg{VcOptim}()(false))
        end
    end
end

end
