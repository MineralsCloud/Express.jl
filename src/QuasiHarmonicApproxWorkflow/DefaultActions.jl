module DefaultActions

using AbInitioSoftwareBase.Inputs: Input, writetxt
using PyQHA: converter, runcode, plot
using SimpleWorkflow: InternalAtomicJob

using ...Express: Action, loadconfig
using ..QuasiHarmonicApproxWorkflow: QuasiHarmonicApprox
using ..Config: materialize
import ..QuasiHarmonicApproxWorkflow: buildjob

struct MakeInput{T} <: Action{T} end
function (x::MakeInput{QuasiHarmonicApprox})(inp_file_list, inp_static, inp_q_points)
    converter(inp_file_list, inp_static, inp_q_points)
end
function (x::MakeInput{QuasiHarmonicApprox})(cfgfile)
    input, file, inp_file_list, static, q_points = loadconfig(cfgfile)
    cd(dirname(input)) do
        x(inp_file_list, static, q_points)
    end
end

buildjob(x::MakeInput, cfgfile) = InternalAtomicJob(() -> x(cfgfile))

struct CalculateThermodyn{T} <: Action{T} end
function (x::CalculateThermodyn{QuasiHarmonicApprox})(cfgfile)
    config = loadconfig(cfgfile)
    runcode(config[:config])
end

buildjob(x::CalculateThermodyn, cfgfile) = InternalAtomicJob(() -> x(cfgfile))

struct Plot{T} <: Action{T} end
function (x::Plot{QuasiHarmonicApprox})(cfgfile)
    config = loadconfig(cfgfile)
    plot(config[:config])
end

buildjob(x::Plot, cfgfile) = InternalAtomicJob(() -> x(cfgfile))

end
