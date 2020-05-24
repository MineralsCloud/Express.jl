using Distributed
addprocs(8)
@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere using Express, Express.EosFitting, Express.Jobs, Express.CLI
@everywhere using EquationsOfState.NonlinearFitting,
    EquationsOfState.Collections, EquationsOfState.Find
@everywhere using QuantumESPRESSO.Inputs.PWscf, QuantumESPRESSOBase.CLI
@everywhere using QuantumESPRESSOParsers
@everywhere using Unitful, UnitfulAtomic
@everywhere using DockerPy.Client, DockerPy.Images, DockerPy.Containers

@everywhere str = raw"""
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

@everywhere docker = DockerClient()
@everywhere image = pull(docker, "rinnocente/qe-full-6.2.1")[1]
# @everywhere begin
#     container = Container(
#         docker,
#         image,
#         command = "sh",
#         name = "qe",
#         tty = true,
#         stdin_open = true,
#         volumes = Dict(
#             joinpath(pwd(), "examples") =>
#                 Dict("bind" => "/home/qe/test", "mode" => "rw"),
#         ),
#     )
# end
@everywhere container = containers(docker)[1]
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
@everywhere pressures = [-50, 0, 50, 100, 150, 200, 250, 300] * u"kbar"
@everywhere scfdirs_local = map(x -> mkpath("examples/scf$(ustrip(x))"), pressures)
@everywhere scfinputs_local = map(x -> x * "/scf.in", scfdirs_local)
@everywhere crude_eos = BirchMurnaghan3rd(317 * u"bohr^3", 210 * u"GPa", 4, -612.43 * u"Ry")
@everywhere template = parse(PWInput, str)
Step(1)(scfinputs_local, template, crude_eos, pressures)
@everywhere scfinputs_docker = map(x -> "/home/qe/test/scf$(ustrip(x))/scf.in", pressures)
@everywhere scfoutputs_docker = map(x -> replace(x, ".in" => ".out"), scfinputs_docker)
# ================================================================= Step 2 =============================
@everywhere bag = Step(2)(
    scfinputs_docker,
    scfoutputs_docker,
    16,
    pwcmd(bin = "/home/qe/qe-6.2.1/bin/pw.x"),
    workers(),
    isdocker = true,
    container = container,
)
map(bag, map(x -> replace(x, ".in" => ".out"), scfinputs_local)) do x, outfile
    touch(outfile)
    open(outfile, "r+") do io
        write(io, x.output[1])
    end
end
# ================================================================= Step 3: read scf.out and curve-fitting =============================
new_eos = Step(3)(map(x -> replace(x, ".in" => ".out"), scfinputs_local), crude_eos)
# new eos:  317.75905077576425 a₀^3, 172.89506496025282 GPa, 4.357510886414555, -612.4315102681139 Ry
# ================================================================= Step 4 =============================
vcdirs_local = map(x -> mkpath("examples/vc$(ustrip(x))"), pressures)
vcinputs_local = map(x -> x * "/vc.in", vcdirs_local)
Step(4)(vcinputs_local, template, new_eos, pressures)
# ================================================================= Step 5 =============================
vcinputs_docker = map(x -> replace(x, "scf" => "vc"), scfinputs_docker)
vcoutput_docker = map(x -> replace(x, ".in" => ".out"), vcinputs_local)
bag2 = Step(5)(
    vcinputs_docker,
    vcoutput_docker,
    16,
    pwcmd(bin = "/home/qe/qe-6.2.1/bin/pw.x"),
    workers(),
    isdocker = true,
    container = container,
)
map(bag2, map(x -> replace(x, ".in" => ".out"), vcinputs_local)) do x, outfile
    touch(outfile)
    open(outfile, "r+") do io
        write(io, x.output[1])
    end
end
stop(container)
# ================================================================= Step 6: read vcrelax.out file and curve-fitting =============================
# final eos: 317.7711705399742 a₀^3, 172.8125730578396 GPa, 4.3535165337769195, -612.4315134858152 Ry
result = Step(6)(vcoutput_docker, new_eos)
println(result)
