struct MakeInput{T} <: Action{T} end
function (::MakeInput{T})(template::S, args...)::S where {T,S<:Input}
    return adjust(template, T(), args...)
end
function (x::MakeInput)(file, template::Input, args...)
    input = x(template, args...)
    mkpath(dirname(file))  # In case its parent directory is not created
    writetxt(file, input)
    return input
end
function (x::MakeInput{T})(cfgfile; kwargs...) where {T}
    config = loadconfig(cfgfile)
    infiles = first.(iofiles(T(), cfgfile))
    eos = PressureEquation(T <: Scf ? config.trial_eos : FitEos{Scf}()(cfgfile))
    return broadcast(
        x,
        infiles,
        config.templates,
        config.pressures,
        fill(eos, length(infiles));
        kwargs...,
    )
end

buildjob(x::MakeInput, cfgfile) = InternalAtomicJob(() -> x(cfgfile))

function adjust end
