using AbInitioSoftwareBase: load
using PyQHA: converter, runcode, plot
using SimpleWorkflows: AtomicJob

using ...Express: Action
using ..QuasiHarmonicApproxWorkflow: QuasiHarmonicApprox
using .Config: ExpandConfig

struct MakeInput{T} <: Action{T} end
function (x::MakeInput{QuasiHarmonicApprox})(inp_file_list, inp_static, inp_q_points)
    converter(inp_file_list, inp_static, inp_q_points)
end

function buildjob(x::MakeInput{QuasiHarmonicApprox}, cfgfile)
    dict = load(cfgfile)
    config = ExpandConfig{QuasiHarmonicApprox}()(dict)
    return AtomicJob(function ()
        return cd(dirname(config.input)) do
            x(config.inp_file_list, config.static, config.q_points)
        end
    end)
end

struct CalculateThermodyn{T} <: Action{T} end

function buildjob(::CalculateThermodyn{QuasiHarmonicApprox}, cfgfile)
    dict = load(cfgfile)
    config = ExpandConfig{QuasiHarmonicApprox}()(dict)
    return AtomicJob(() -> runcode(config.config))
end

struct Plot{T} <: Action{T} end

function buildjob(::Plot{QuasiHarmonicApprox}, cfgfile)
    dict = load(cfgfile)
    config = ExpandConfig{QuasiHarmonicApprox}()(dict)
    return AtomicJob(() -> plot(config.config))
end
