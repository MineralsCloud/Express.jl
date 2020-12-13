struct MakeInput{T} <: Action{T} end
function (::MakeInput{T})(template::S, args...)::S where {T,S<:Input}
    return adjust(template, T(), args...)
end
function (x::MakeInput)(file, template::Input, args...)
    input = x(template, args...)
    mkpath(dirname(file))  # In case its parent directory is not created
    writeinput(file, input)
    return input
end
function (x::MakeInput{T})(cfgfile; kwargs...) where {T}
    settings = loadconfig(cfgfile)
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

buildjob(::MakeInput{T}, cfgfile) where {T} =
    InternalAtomicJob(() -> MakeInput{T}()(cfgfile))

function adjust end
