struct FitEos{T} <: Action{T} end
function (x::FitEos{T})(
    data::AbstractVector{<:Pair},
    trial_eos::EnergyEOS,
) where {T<:ScfOrOptim}
    return eosfit(trial_eos, first.(data), last.(data))
end
function (x::FitEos{T})(outputs, trial_eos::EnergyEOS) where {T<:ScfOrOptim}
    data = GetData{T}()(outputs)
    if length(data) <= 5
        @info "pressures <= 5 may give unreliable results, run more if possible!"
    end
    return x(data, trial_eos)
end
function (x::FitEos{T})(cfgfile) where {T<:ScfOrOptim}
    settings = load_settings(cfgfile)
    outfiles = last.(iofiles(T(), cfgfile))
    rawsettings = load(cfgfile)
    saveto = joinpath(rawsettings["workdir"], shortname(T) * "_eos.jls")
    trial_eos = T <: Scf ? settings.trial_eos : deserialize(saveto)
    eos = x(outfiles, EnergyEOS(trial_eos))
    SaveEos{T}(saveto, eos)
    return eos
end
