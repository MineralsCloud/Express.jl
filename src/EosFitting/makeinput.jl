struct MakeInput{T} <: Action{T}
    calc::T
end
function (x::MakeInput)(file, template::Input, args...)
    input = customize(standardize(template, x.calc), args...)
    writeinput(file, input)
    return input
end
function (x::MakeInput)(cfgfile; kwargs...)
    settings = load_settings(cfgfile)
    infiles = first.(iofiles(x.calc, cfgfile))
    eos = PressureEOS(x.calc isa Scf ? settings.trial_eos : FitEos(Scf())(cfgfile))
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
