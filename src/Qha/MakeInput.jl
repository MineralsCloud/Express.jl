struct MakeInput{T} <: Action{T} end
function (x::MakeInput{QuasiHarmonicApprox})(inp_file_list, inp_static, inp_q_points)
    converter(inp_file_list, inp_static, inp_q_points)
end
function (x::MakeInput{QuasiHarmonicApprox})(cfgfile)
    input, file, inp_file_list, static, q_points = loadconfig(cfgfile)
    x(inp_file_list, static, q_points)
end

buildjob(x::MakeInput, cfgfile) = InternalAtomicJob(() -> x(cfgfile))
