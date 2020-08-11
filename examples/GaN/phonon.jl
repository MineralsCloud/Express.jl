using Distributed
addprocs(16)
using Express.Phonon, Express.Jobs
using QuantumESPRESSO.Inputs.PWscf, QuantumESPRESSO.Inputs.PHonon, QuantumESPRESSO.CLI
using QuantumESPRESSOParsers
using Unitful, UnitfulAtomic

pressures = [0, 50] * u"kbar"
vcoutdirs_local = map(x -> "/Users/qz/Downloads/test/vc$(ustrip(x))/", pressures)
vcout_local = map(x -> x * "vc.out", vcoutdirs_local)
template = parse(PWInput, read("examples/GaN/template.in", String))
scfdirs = map(x -> "examples/GaN/scf$x", pressures)
scfinputs = map(x -> x * "/scf.in", scfdirs)
preprocess(DfptMethod(), scfinputs[1], vcout_local[1], template)
# ================================================================= Step 2 =============================
bag = Step(2)(
    scfinputs_docker,
    scfoutputs_docker,
    4,
    pwcmd(bin = "/home/qe/qe-6.2.1/bin/pw.x"),
    [2, 3],
    isdocker = true,
    container = container,
)
phononinputs_local = map(x -> x * "/ph.in", scfdirs_local)
preprocess(
    Step(2),
    phononinputs_local,
    scfinputs_local,
    PhInput(PhNamelist(nq1 = 2, nq2 = 2, nq3 = 2, tr2_ph = 1e-14, ldisp = true)),
)
bag2 = Step(2)(
    phononinputs_docker,
    phononoutputs_docker,
    16,
    pwcmd(bin = "/home/qe/qe-6.2.1/bin/ph.x"),
    [2, 3],
    isdocker = true,
    container = container,
)
map(bag2, map(x -> replace(x, ".in" => ".out"), phononinputs_local)) do x, outfile
    touch(outfile)
    open(outfile, "w") do io
        write(io, join(p[1] for p in [s for s in x.output]))
    end
end
