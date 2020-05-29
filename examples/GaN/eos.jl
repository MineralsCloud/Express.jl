using Distributed
addprocs(6)
using Pkg
Pkg.activate(".")
using Express, Express.EosFitting, Express.Jobs, Express.CLI, Express.Workspaces
using EquationsOfState.NonlinearFitting, EquationsOfState.Collections, EquationsOfState.Find
using QuantumESPRESSO.Inputs.PWscf, QuantumESPRESSOBase.CLI
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
container = Container(
    docker,
    image,
    command = "sh",
    name = "qe",
    tty = true,
    stdin_open = true,
    volumes = Dict(
        joinpath(pwd(), "examples/GaN") => Dict("bind" => "/home/qe/test", "mode" => "rw"),
    ),
)
# container = containers(docker)[1]
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
# ================================================================= Step 1 =============================
# pressures = [-100, -50, 0, 50, 100, 150, 200, 250, 300, 350, 400, 500] * u"kbar"
pressures = [-50, 0, 50, 100, 200, 300] * u"kbar"
scfdirs_local = map(x -> mkpath("examples/GaN/scf$(ustrip(x))"), pressures)
scfinputs_local = map(x -> x * "/scf.in", scfdirs_local)
crude_eos = BirchMurnaghan3rd(317 * u"bohr^3", 210 * u"GPa", 4, -612.43 * u"Ry")
template = parse(PWInput, str)
SelfConsistentField{PREPARE_INPUT}()(scfinputs_local, template, crude_eos, pressures)
scfinputs_docker = map(x -> "/home/qe/test/scf$(ustrip(x))/scf.in", pressures)
scfoutputs = map(x -> replace(x, ".in" => ".out"), scfinputs_local)
# ================================================================= Step 2 =============================
bag = SelfConsistentField{LAUNCH_JOB}()(
    scfinputs_docker,
    scfoutputs,
    DockerWorkspace(12, container),
    pwcmd(bin = "/home/qe/qe-6.2.1/bin/pw.x"),
    workers(),
)
# ================================================================= Step 3: read scf.out and curve-fitting =============================
new_eos = SelfConsistentField{ANALYSE_OUTPUT}()(
    map(x -> replace(x, ".in" => ".out"), scfinputs_local),
    crude_eos,
)
# new eos:  317.75905077576425 a₀^3, 172.89506496025282 GPa, 4.357510886414555, -612.4315102681139 Ry
# ================================================================= Step 4 =============================
vcdirs_local = map(x -> mkpath("examples/GaN/vc$(ustrip(x))"), pressures)
vcinputs_local = map(x -> x * "/vc.in", vcdirs_local)
SelfConsistentField{PREPARE_INPUT}()(vcinputs_local, template, new_eos, pressures)
# ================================================================= Step 5 =============================
vcinputs_docker = map(x -> replace(x, "scf" => "vc"), scfinputs_docker)
vcoutput = map(x -> replace(x, ".in" => ".out"), vcinputs_local)
bag2 = SelfConsistentField{LAUNCH_JOB}()(
    vcinputs_docker,
    vcoutput,
    DockerWorkspace(12, container),
    pwcmd(bin = "/home/qe/qe-6.2.1/bin/pw.x"),
    workers(),
)
stop(container)
# ================================================================= Step 6: read vcrelax.out file and curve-fitting =============================
# final eos: 317.7711705399742 a₀^3, 172.8125730578396 GPa, 4.3535165337769195, -612.4315134858152 Ry
result = SelfConsistentField{ANALYSE_OUTPUT}()(vcoutput, new_eos)
println(result)
