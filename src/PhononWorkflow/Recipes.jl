module Recipes

using AbInitioSoftwareBase: load
using ExpressBase: Scf, Dfpt, RealSpaceForceConstants, PhononDispersion, VDos
using ExpressBase.Recipes: Recipe
using SimpleWorkflows.Jobs: Job
using SimpleWorkflows.Workflows: Workflow, run!, →, ⇉, ⇶, ⭃

using ...Express: DownloadPotentials, LogTime, RunCmd, jobify
using ..PhononWorkflow: MakeInput

struct PhononDispersionRecipe <: Recipe
    config
end
struct VDosRecipe <: Recipe
    config
end

function build(::Type{Workflow}, r::PhononDispersionRecipe)
    a0 = jobify(DownloadPotentials{Scf}(), r.config)
    a = jobify(LogTime{Scf}())
    b = jobify(MakeInput{Scf}(), r.config)
    c = jobify(RunCmd{Scf}(), r.config)
    f = jobify(MakeInput{Dfpt}(), r.config)
    g = jobify(RunCmd{Dfpt}(), r.config)
    i = jobify(LogTime{RealSpaceForceConstants}())
    j = jobify(MakeInput{RealSpaceForceConstants}(), r.config)
    k = jobify(RunCmd{RealSpaceForceConstants}(), r.config)
    m = jobify(LogTime{PhononDispersion}())
    n = jobify(MakeInput{PhononDispersion}(), r.config)
    o = jobify(RunCmd{PhononDispersion}(), r.config)
    a0 → a ⇉ b ⇶ c ⇉ f ⇶ g ⭃ i ⇉ j ⇶ k ⭃ m ⇉ n ⇶ o
    return Workflow(a0)
end
function build(::Type{Workflow}, r::VDosRecipe)
    a0 = jobify(DownloadPotentials{Scf}(), r.config)
    a = jobify(LogTime{Scf}())
    b = jobify(MakeInput{Scf}(), r.config)
    c = jobify(RunCmd{Scf}(), r.config)
    f = jobify(MakeInput{Dfpt}(), r.config)
    g = jobify(RunCmd{Dfpt}(), r.config)
    i = jobify(LogTime{RealSpaceForceConstants}())
    j = jobify(MakeInput{RealSpaceForceConstants}(), r.config)
    k = jobify(RunCmd{RealSpaceForceConstants}(), r.config)
    m = jobify(LogTime{VDos}())
    n = jobify(MakeInput{VDos}(), r.config)
    o = jobify(RunCmd{VDos}(), r.config)
    a0 → a ⇉ b ⇶ c ⇉ f ⇶ g ⭃ i ⇉ j ⇶ k ⭃ m ⇉ n ⇶ o
    return Workflow(a0)
end
function build(::Type{Workflow}, file)
    dict = load(file)
    recipe = dict["recipe"]
    if recipe == "phonon dispersion"
        return build(Workflow, PhononDispersionRecipe(dict))
    elseif dict["recipe"] == "vdos"
        return build(Workflow, VDosRecipe(dict))
    else
        error("unsupported recipe $recipe.")
    end
end

end
