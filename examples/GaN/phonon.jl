using Distributed
addprocs(16)
using Express.Phonon, Express.Jobs, Express.CLI
using QuantumESPRESSO.Inputs.PWscf, QuantumESPRESSO.Inputs.PHonon, QuantumESPRESSOBase.CLI
using QuantumESPRESSOParsers
using Unitful, UnitfulAtomic
using DockerPy.Client, DockerPy.Images, DockerPy.Containers

str = raw"""
&control
calculation='scf'
verbosity = 'high'
pseudo_dir  = '/home/qe/pseudo'
prefix='GaN'
outdir      = './tmp'
tprnfor = .true.
tstress = .true.
/
&system
ibrav = 4
celldm(1)=5.95484286816
celldm(3)=1.63011343669
nat = 4
ntyp = 2
ecutwfc = 160
/
&electrons
diagonalization='david'
mixing_beta=0.7
conv_thr=1.0d-10
/
ATOMIC_SPECIES
Ga  0   Ga.pbe-dn-kjpaw_psl.1.0.0.UPF
N   0   N.pbe-n-kjpaw_psl.1.0.0.UPF
ATOMIC_POSITIONS (crystal)
Ga       0.666666667   0.333333333  -0.000051966
N        0.666666667   0.333333333   0.376481188
Ga       0.333333333   0.666666667   0.499948034
N        0.333333333   0.666666667   0.876481188
K_POINTS automatic
6 6 6 0 0 0
"""

docker = DockerClient()
image = pull(docker, "rinnocente/qe-full-6.2.1")[1]
try
    container = Container(
        docker,
        image,
        command = "sh",
        name = "qe",
        tty = true,
        stdin_open = true,
        volumes = Dict(
            joinpath(pwd(), "examples") =>
                Dict("bind" => "/home/qe/test", "mode" => "rw"),
        ),
    )
catch
    container = get(docker.containers, "qe")
end
start(container)
exec_run(container, "mkdir -p /home/qe/pseudo/")
exec_run(
    container,
    "wget -O /home/qe/pseudo/Ga.pbe-dn-kjpaw_psl.1.0.0.UPF https://www.quantum-espresso.org/upf_files/Ga.pbe-dn-kjpaw_psl.1.0.0.UPF",
)
exec_run(
    container,
    "wget -O /home/qe/pseudo/N.pbe-n-kjpaw_psl.1.0.0.UPF https://www.quantum-espresso.org/upf_files/N.pbe-n-kjpaw_psl.1.0.0.UPF",
)

pressures = [0, 50] * u"kbar"
vcoutdirs_local = map(x -> "/Users/qz/Downloads/test/vc$(ustrip(x))/", pressures)
vcout_local = map(x -> x * "vc.out", vcoutdirs_local)
scfdirs_local = map(x -> mkpath("examples/scf$(ustrip(x))"), pressures)
scfinputs_local = map(x -> x * "/scf.in", scfdirs_local)
template = parse(PWInput, str)
Step(1)(scfinputs_local, vcout_local, template)
scfinputs_docker = map(x -> "/home/qe/test/scf$(ustrip(x))/scf.in", pressures)
scfoutputs_docker = map(x -> replace(x, ".in" => ".out"), scfinputs_docker)
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
map(bag, map(x -> replace(x, ".in" => ".out"), scfinputs_local)) do x, outfile
    touch(outfile)
    open(outfile, "w") do io
        write(io, join(p[1] for p in [s for s in x.output]))
    end
end
# ===
phononinputs_local = map(x -> x * "/ph.in", scfdirs_local)
phononinputs_docker = map(x -> "/home/qe/test/scf$(ustrip(x))/ph.in", pressures)
phononoutputs_docker = map(x -> replace(x, ".in" => ".out"), phononinputs_docker)
preprocess(Step(2), phononinputs_local, scfinputs_local, PhInput(PhNamelist(nq1 = 2, nq2 = 2, nq3 = 2, tr2_ph = 1e-14, ldisp = true)))
# ===
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
