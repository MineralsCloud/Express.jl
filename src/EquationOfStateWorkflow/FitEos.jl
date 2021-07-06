struct FitEos{T} <: Action{T} end
function (x::FitEos{T})(
    data::AbstractVector{<:Pair},
    trial_eos::EnergyEquation,
) where {T<:ScfOrOptim}
    return eosfit(trial_eos, first.(data), last.(data))
end
function (x::FitEos{T})(outputs, trial_eos::EnergyEquation) where {T<:ScfOrOptim}
    data = GetData{T}()(outputs)
    if length(data) <= 5
        @info "pressures <= 5 may give unreliable results, run more if possible!"
    end
    return x(data, trial_eos)
end
function (x::FitEos{T})(cfgfile) where {T<:ScfOrOptim}
    config = loadconfig(cfgfile)
    outfiles = last.(iofiles(T(), cfgfile))
    saveto = joinpath(config.workdir, shortname(T) * "_eos.jls")
    trial_eos =
        T <: Scf ? config.trial_eos :
        deserialize(joinpath(config.workdir, shortname(Scf) * "_eos.jls"))
    eos = x(outfiles, EnergyEquation(trial_eos))
    SaveEos{T}()(saveto, eos)
    return eos
end
