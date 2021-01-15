struct CalculateThermodyn{T} <: Action{T} end
function (x::CalculateThermodyn{QuasiHarmonicApprox})(cfgfile)
    config = loadconfig(cfgfile)
    runcode(config[:config])
end
