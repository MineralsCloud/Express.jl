struct MakeInput{T} <: Action{T} end
function (x::MakeInput{T})(file, template::S, args...)::S where {T,S<:Input}
    update = customize(args...) ∘ standardize(T())
    input = update(template)
    writeinput(file, input)
    return input
end
function (x::MakeInput{T})(cfgfile; kwargs...) where {T}
    settings = load_settings(cfgfile)
    infiles = first.(iofiles(T(), cfgfile))
    eos = PressureEOS(T <: Scf ? settings.trial_eos : FitEos{Scf}()(cfgfile))
    return broadcast(
        x,
        infiles,
        settings.templates,
        settings.pressures_or_volumes,
        fill(eos, length(infiles));
        kwargs...,
    )
end

function standardize end

function customize end
