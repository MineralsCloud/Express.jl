struct Plot{T} <: Action{T} end
function (x::Plot{QuasiHarmonicApprox})(cfgfile)
    config = loadconfig(cfgfile)
    plot(config[:config])
end
