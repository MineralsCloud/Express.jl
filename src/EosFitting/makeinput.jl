struct MakeInput{T} <: Action{T} end
MakeInput(::T) where {T<:Calculation} = MakeInput{T}()
function (::MakeInput{T})(file, template::Input, args...) where {T<:ScfOrOptim}
    input = customize(standardize(template, T()), args...)
    writeinput(file, input)
    return input
end
function (x::MakeInput{T})(cfgfile; kwargs...) where {T<:ScfOrOptim}
    settings = load_settings(cfgfile)
    infiles = first.(iofiles(T(), cfgfile))
    eos = PressureEOS(
        T == SelfConsistentField ? settings.trial_eos :
        FitEos(SelfConsistentField())(cfgfile),
    )
    return broadcast(
        x,
        infiles,
        settings.templates,
        settings.pressures_or_volumes,
        fill(eos, length(infiles));
        kwargs...,
    )
end
const makeinput = MakeInput

function standardize end

function customize end
