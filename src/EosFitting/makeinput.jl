struct MakeInput{T} <: Action{T} end
function (x::MakeInput{T})(file, template::Input, args...) where {T}
    modify = customize(args...) ∘ standardize(T())
    input = modify(template)
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
