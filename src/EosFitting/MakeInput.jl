struct MakeInput{T} <: Action{T} end
function (x::MakeInput)(file, template::Input, args...)
    input = x(template, args...)
    dir = dirname(file)
    if !isdir(dir)
        @warn "directory `$dir` not found! I will create it for you."
        mkpath(dir)
    end
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
