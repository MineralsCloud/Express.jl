using AbInitioSoftwareBase: load
using ExpressBase: QuasiHarmonicApproximation, Action
using PyQHA: converter, runcode, plot
using SimpleWorkflows.Jobs: Job

using .Config: ExpandConfig

import ExpressWorkflowMaker.Templates: jobify

struct MakeInput{T} <: Action{T} end
function (x::MakeInput{QuasiHarmonicApproximation})(inp_file_list, inp_static, inp_q_points)
    converter(inp_file_list, inp_static, inp_q_points)
end

function jobify(x::MakeInput{QuasiHarmonicApproximation}, cfgfile)
    dict = load(cfgfile)
    config = ExpandConfig{QuasiHarmonicApproximation}()(dict)
    return Job(function ()
        return cd(dirname(config.input)) do
            x(config.inp_file_list, config.static, config.q_points)
        end
    end)
end

struct CalculateThermodyn{T} <: Action{T} end

function jobify(::CalculateThermodyn{QuasiHarmonicApproximation}, cfgfile)
    dict = load(cfgfile)
    config = ExpandConfig{QuasiHarmonicApproximation}()(dict)
    return Job(() -> runcode(config.config))
end

struct Plot{T} <: Action{T} end

function jobify(::Plot{QuasiHarmonicApproximation}, cfgfile)
    dict = load(cfgfile)
    config = ExpandConfig{QuasiHarmonicApproximation}()(dict)
    return Job(() -> plot(config.config))
end
