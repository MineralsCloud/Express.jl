module DefaultActions

using AbInitioSoftwareBase: load
using AbInitioSoftwareBase.Inputs: Input, writetxt
using Dates: now, format
using Logging: with_logger, current_logger
using SimpleWorkflows: AtomicJob
using ...Express: Action, Calculation, LatticeDynamics, Scf, calculation
using ...Shell: distprocs
using ...EquationOfStateWorkflow: VcOptim
using ..PhononWorkflow: Dfpt, RealSpaceForceConstants, PhononDispersion, VDos
using ..Config: ExpandConfig

struct MakeInput{T} <: Action{T} end
function (x::MakeInput{T})(
    file,
    template::Input,
    args...,
) where {T<:Union{Scf,LatticeDynamics}}
    input = x(template, args...)
    mkpath(dirname(file))  # In case its parent directory is not created
    writetxt(file, input)
    return input
end

function buildjob(x::MakeInput{Dfpt}, cfgfile)
    dict = load(cfgfile)
    pop!(dict, "workflow")
    config = ExpandConfig{Dfpt}()(dict)
    return map(config.files[2], config.files[1]) do (input, _), (previnput, _)
        AtomicJob(() -> x(input, config.template[2], previnput))
    end
end
function buildjob(x::MakeInput{RealSpaceForceConstants}, cfgfile)
    dict = load(cfgfile)
    pop!(dict, "workflow")
    config = ExpandConfig{RealSpaceForceConstants}()(dict)
    return map(config.files[3], config.files[2]) do (input, _), (previnput, _)
        AtomicJob(() -> x(input, config.template[3], previnput))
    end
end
function buildjob(x::MakeInput{Scf}, cfgfile)
    dict = load(cfgfile)
    pop!(dict, "workflow")
    config = ExpandConfig{Scf}()(dict)
    cells = map(config.files[1]) do (input, _)
        file = joinpath(dirname(input), string(VcOptim()) * ".out")
        str = read(file, String)
        cell = parsecell(str)
        if any(x === nothing for x in cell)
            error("set cell failed!")
        else
            cell
        end
    end
    return map(config.files[1], cells) do file, cell
        input = first(file)
        AtomicJob(() -> x(input, config.template[1], first(cell), last(cell)))
    end
end
function buildjob(x::MakeInput{T}, cfgfile) where {T<:Union{PhononDispersion,VDos}}
    dict = load(cfgfile)
    pop!(dict, "workflow")
    config = ExpandConfig{T}()(dict)
    ifcinputs = map(config.files[3]) do (input, _)
        parse(inputtype(RealSpaceForceConstants), read(input, String))
    end
    dfptinputs = map(config.files[2]) do (input, _)
        parse(inputtype(Dfpt), read(input, String))
    end
    return map(config.files[4], ifcinputs, dfptinputs) do (input, _), ifcinput, dfptinput
        AtomicJob(() -> x(input, config.template[2], ifcinput, dfptinput))
    end
end

function parsecell end

function inputtype end

struct RunCmd{T} <: Action{T} end

function buildjob(x::RunCmd{Scf}, cfgfile)
    dict = load(cfgfile)
    pop!(dict, "workflow")
    config = ExpandConfig{Scf}()(dict)
    np = distprocs(config.cli.mpi.np, length(config.files[1]))
    return map(config.files[1]) do (input, output)
        AtomicJob(() -> x(input, output; np = np))
    end
end
function buildjob(x::RunCmd{Dfpt}, cfgfile)
    dict = load(cfgfile)
    pop!(dict, "workflow")
    config = ExpandConfig{Dfpt}()(dict)
    np = distprocs(config.cli.mpi.np, length(config.files[2]))
    return map(config.files[2]) do (input, output)
        AtomicJob(() -> x(input, output; np = np))
    end
end
function buildjob(x::RunCmd{RealSpaceForceConstants}, cfgfile)
    dict = load(cfgfile)
    pop!(dict, "workflow")
    config = ExpandConfig{RealSpaceForceConstants}()(dict)
    np = distprocs(config.cli.mpi.np, length(config.files[3]))
    return map(config.files[3]) do (input, output)
        AtomicJob(() -> x(input, output; np = np))
    end
end
function buildjob(x::RunCmd{T}, cfgfile) where {T<:LatticeDynamics}
    dict = load(cfgfile)
    pop!(dict, "workflow")
    config = ExpandConfig{T}()(dict)
    np = distprocs(config.cli.mpi.np, length(config.files[4]))
    return map(config.files[4]) do (input, output)
        AtomicJob(() -> x(input, output; np = np))
    end
end

struct LogMsg{T} <: Action{T} end
function (x::LogMsg)(; start = true)
    act = start ? "starts" : "ends"
    with_logger(current_logger()) do
        println(
            "The calculation $(calculation(x)) $act at $(format(now(), "HH:MM:SS u dd, yyyy")).",
        )
    end
end

end
