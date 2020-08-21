using EquationsOfStateOfSolids.Collections
using Express.EosFitting
using QuantumESPRESSO.Inputs.PWscf
using Unitful, UnitfulAtomic

# pressures = [-100, -50, 0, 50, 100, 150, 200, 250, 300, 350, 400, 500] * u"kbar"
pressures = [-50, 0, 50, 100, 200, 300] * u"kbar"
crude_eos = BirchMurnaghan3rd(317 * u"bohr^3", 210 * u"GPa", 4, -612.43 * u"Ry")
template = parse(PWInput, read("examples/GaN/template.in", String))
# ================================================================= Step 1 =============================
scfdirs = map(x -> "examples/GaN/scf$x", pressures)
scfinputs = map(x -> x * "/scf.in", scfdirs)
SelfConsistentField()(PREPARE_INPUT)(scfinputs, template, pressures, crude_eos)
# ================================================================= Step 2 =============================
scfoutputs = map(x -> replace(x, ".in" => ".out"), scfinputs)
bag = SelfConsistentField()(LAUNCH_JOB)(
    scfoutputs,
    scfinputs,
    12,
    "/home/qe/qe-6.2.1/bin/pw.x",
)
# ================================================================= Step 3: read scf.out and curve-fitting =============================
new_eos = SelfConsistentField()(ANALYSE_OUTPUT)(scfoutputs, crude_eos)
# new eos:  317.75905077576425 a₀^3, 172.89506496025282 GPa, 4.357510886414555, -612.4315102681139 Ry
# ================================================================= Step 4 =============================
vcdirs = map(x -> mkpath("examples/GaN/vc$x"), pressures)
vcinputs = map(x -> x * "/vc.in", vcdirs)
VariableCellOptimization()(PREPARE_INPUT)(vcinputs, template, pressures, new_eos)
# ================================================================= Step 5 =============================
vcoutputs = map(x -> replace(x, ".in" => ".out"), vcinputs)
bag2 = VariableCellOptimization()(LAUNCH_JOB)(
    vcoutputs,
    vcinputs,
    12,
    "/home/qe/qe-6.2.1/bin/pw.x",
)
# ================================================================= Step 6: read vcrelax.out file and curve-fitting =============================
# final eos: 317.7711705399742 a₀^3, 172.8125730578396 GPa, 4.3535165337769195, -612.4315134858152 Ry
result = VariableCellOptimization()(ANALYSE_OUTPUT)(vcoutputs, new_eos)
println(result)
