using AbInitioSoftwareBase: load
using AbInitioSoftwareBase.Inputs: Input, writetxt, getpseudodir, getpotentials
using Dates: now, format
using ExpressBase:
    Calculation,
    LatticeDynamics,
    Scf,
    Dfpt,
    RealSpaceForceConstants,
    PhononDispersion,
    VDos,
    Action,
    calculation
using SimpleWorkflows: Job

using ..Shell: distprocs
using ..EquationOfStateWorkflow: VariableCellOptimization
using .Config: ExpandConfig

import ExpressWorkflowMaker.Templates: jobify

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

function jobify(x::MakeInput{Dfpt}, cfgfile)
    dict = load(cfgfile)
    config = ExpandConfig{Dfpt}()(dict)
    return map(config.files[2], config.files[1]) do (input, _), (previnput, _)
        Job(function ()
            str = read(previnput, String)
            prev = parse(inputtype(Scf), str)
            return x(input, config.template.dfpt, prev)
        end)
    end
end
function jobify(x::MakeInput{RealSpaceForceConstants}, cfgfile)
    dict = load(cfgfile)
    config = ExpandConfig{RealSpaceForceConstants}()(dict)
    return map(config.files[3], config.files[2]) do (input, _), (previnput, _)
        Job(function ()
            str = read(previnput, String)
            prev = parse(inputtype(Dfpt), str)
            return x(input, config.template.q2r, prev)
        end)
    end
end
function jobify(x::MakeInput{T}, cfgfile) where {T<:Union{PhononDispersion,VDos}}
    dict = load(cfgfile)
    config = ExpandConfig{T}()(dict)
    return map(
        config.files[4],
        config.files[3],
        config.files[2],
    ) do (input, _), (previnput, _), (pprevinput, _)
        Job(function ()
            str = read(previnput, String)
            q2r = parse(inputtype(RealSpaceForceConstants), str)
            str = read(pprevinput, String)
            dfpt = parse(inputtype(Dfpt), str)
            return x(input, config.template.disp, q2r, dfpt)
        end)
    end
end
function jobify(x::MakeInput{Scf}, cfgfile)
    dict = load(cfgfile)
    config = ExpandConfig{Scf}()(dict)
    return map(config.files[1]) do (input, _)
        Job(
            function ()
                vcout = joinpath(
                    dirname(input),
                    string(nameof(VariableCellOptimization)) * ".out",
                )
                str = read(vcout, String)
                cell = parsecell(str)
                if any(x === nothing for x in cell)
                    error("set cell failed!")
                else
                    cell
                end
                return x(input, config.template.scf, first(cell), last(cell))
            end,
        )
    end
end

function parsecell end

function inputtype end

function jobify(x::RunCmd{Scf}, cfgfile)
    dict = load(cfgfile)
    config = ExpandConfig{Scf}()(dict)
    np = distprocs(config.cli.mpi.np, length(config.files[1]))
    return map(config.files[1]) do (input, output)
        Job(() -> x(input, output; np = np))
    end
end
function jobify(x::RunCmd{Dfpt}, cfgfile)
    dict = load(cfgfile)
    config = ExpandConfig{Dfpt}()(dict)
    np = distprocs(config.cli.mpi.np, length(config.files[2]))
    return map(config.files[2]) do (input, output)
        Job(() -> x(input, output; np = np))
    end
end
function jobify(x::RunCmd{RealSpaceForceConstants}, cfgfile)
    dict = load(cfgfile)
    config = ExpandConfig{RealSpaceForceConstants}()(dict)
    np = distprocs(config.cli.mpi.np, length(config.files[3]))
    return map(config.files[3]) do (input, output)
        Job(() -> x(input, output; np = np))
    end
end
function jobify(x::RunCmd{T}, cfgfile) where {T<:LatticeDynamics}
    dict = load(cfgfile)
    config = ExpandConfig{T}()(dict)
    np = distprocs(config.cli.mpi.np, length(config.files[4]))
    return map(config.files[4]) do (input, output)
        Job(() -> x(input, output; np = np))
    end
end
