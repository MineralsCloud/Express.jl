var documenterSearchIndex = {"docs":
[{"location":"examples/GaN/#hcp-GaN-example","page":"hcp-GaN example","title":"hcp-GaN example","text":"","category":"section"},{"location":"examples/GaN/#Fitting-equations-of-state","page":"hcp-GaN example","title":"Fitting equations of state","text":"","category":"section"},{"location":"examples/GaN/#Run-interactively","page":"hcp-GaN example","title":"Run interactively","text":"","category":"section"},{"location":"examples/GaN/","page":"hcp-GaN example","title":"hcp-GaN example","text":"Install the latest version of Julia (as new as you can). Install this package as instructed in \"Installation\" section.","category":"page"},{"location":"examples/GaN/","page":"hcp-GaN example","title":"hcp-GaN example","text":"After the last command, this repo will be cloned to DEPOT_PATH. On a *nix system, it is usually ~/.julia/dev/. Go to this ~/.julia/dev/Express, start your VS Code or Atom at exactly there. Open a Julia REPL in your VS Code/Atom, run","category":"page"},{"location":"examples/GaN/","page":"hcp-GaN example","title":"hcp-GaN example","text":"julia> using Express.EquationOfStateWorkflow.Recipes\n\njulia> config = \"<PATH-TO-EXPRESS>/examples/GaN/eos.yaml\";\n\njulia> buildworkflow(config)","category":"page"},{"location":"examples/GaN/","page":"hcp-GaN example","title":"hcp-GaN example","text":"Then it will start running scf calculations on different pressures. You can type scfjobs in terminal at any time to inquire the status of the jobs. When they are running, a 🚧 emoji will be shown. If succeeded, a ✅ will be shown. If failed, a ❌ will be shown.","category":"page"},{"location":"examples/GaN/","page":"hcp-GaN example","title":"hcp-GaN example","text":"When all jobs are finished (you need at least 6 pressures finished), run","category":"page"},{"location":"examples/GaN/","page":"hcp-GaN example","title":"hcp-GaN example","text":"julia> postprocess(SelfConsistentField(), config)  # Step 3\nEquationsOfState.Collections.BirchMurnaghan3rd{Unitful.Quantity{Float64,D,U} where U where D}\n v0 = 317.7479585715598 a₀³\n b0 = 172.85031797933803 GPa\n b′0 = 4.379649372725796\n e0 = -612.4315064611411 Ry","category":"page"},{"location":"examples/GaN/","page":"hcp-GaN example","title":"hcp-GaN example","text":"An equation of state, with units, will be shown.","category":"page"},{"location":"examples/GaN/","page":"hcp-GaN example","title":"hcp-GaN example","text":"The next 3 steps are basically the same, just with VariableCellOptimization as calculation type.","category":"page"},{"location":"examples/GaN/","page":"hcp-GaN example","title":"hcp-GaN example","text":"julia> preprocess(VariableCellOptimization(), config)  # Step 1\n\njulia> vcjobs = process(VariableCellOptimization(), config)  # Step 2\n\njulia> postprocess(VariableCellOptimization(), config)\nEquationsOfState.Collections.BirchMurnaghan3rd{Unitful.Quantity{Float64,D,U} where U where D}\n v0 = 317.7517433287492 a₀³\n b0 = 172.74897940782353 GPa\n b′0 = 4.388034458575274\n e0 = -612.4315074654779 Ry","category":"page"},{"location":"examples/GaN/","page":"hcp-GaN example","title":"hcp-GaN example","text":"When it is calling Quantum ESPRESSO to do the calculations, it looks like the REPL got stuck, but it is not. Just wait! Do not stop the REPL! Each scf calculation will take about 2 minutes on 2 processors, and each vc-relax calculation will take about 6-9 minutes. So it might need 2-3 hours to run the whole workflow, depending on how good your computer is.","category":"page"},{"location":"examples/GaN/#Run-using-a-script","page":"hcp-GaN example","title":"Run using a script","text":"","category":"section"},{"location":"examples/GaN/","page":"hcp-GaN example","title":"hcp-GaN example","text":"When running on a high performance computer, it is always recommended to first start a debugging/developing environment, do a small test to make sure that inputs can be generated and job can be submitted. Then to scale your calculation, write the above commands in a Julia script and submit it to a schedule manager. That is, write","category":"page"},{"location":"examples/GaN/","page":"hcp-GaN example","title":"hcp-GaN example","text":"using Express.EquationOfStateWorkflow\n\nconfig = \"<PATH-TO-EXPRESS>/examples/GaN/eos.yaml\"\npreprocess(SelfConsistentField(), config)  # Step 1\nscfjobs = process(SelfConsistentField(), config)  # Step 2\nprint(postprocess(SelfConsistentField(), config))  # Step 3\npreprocess(VariableCellOptimization(), config)  # Step 4\nvcjobs = process(VariableCellOptimization(), config)  # Step 5\nprint(postprocess(VariableCellOptimization(), config))  # Step 6","category":"page"},{"location":"examples/GaN/","page":"hcp-GaN example","title":"hcp-GaN example","text":"to a .jl file, say GaN.jl. If you are using Slurm, then write a job.sh with","category":"page"},{"location":"examples/GaN/","page":"hcp-GaN example","title":"hcp-GaN example","text":"#!/bin/sh\n#SBATCH -A <your-account>\n#SBATCH -N 1\n#SBATCH --tasks-per-node=24\n#SBATCH -J job\n#SBATCH --time=02:00:00\n\njulia <PATH-TO-GaN.jl>","category":"page"},{"location":"examples/GaN/","page":"hcp-GaN example","title":"hcp-GaN example","text":"Then submit it using sbatch job.sh.","category":"page"},{"location":"develop/#How-to-develop-this-package","page":"Development","title":"How to develop this package","text":"","category":"section"},{"location":"develop/#Download-the-project","page":"Development","title":"Download the project","text":"","category":"section"},{"location":"develop/","page":"Development","title":"Development","text":"First, add an SSH key to your GitHub account. Request or wait for the administrator of this repo to give you an invitation for triage.","category":"page"},{"location":"develop/","page":"Development","title":"Development","text":"Then, similar to section \"Installation\", run","category":"page"},{"location":"develop/","page":"Development","title":"Development","text":"julia> using Pkg; Pkg.update()\n\njulia> pkg\"dev Express\"","category":"page"},{"location":"develop/","page":"Development","title":"Development","text":"Then the package will be cloned to your local machine at a path. On macOS, by default is located at ~/.julia/dev/Express unless you modify the JULIA_DEPOT_PATH environment variable. (See Julia's official documentation on how to do this.) In the following text, we will call it PKGROOT.","category":"page"},{"location":"develop/#instantiating","page":"Development","title":"Instantiate the project","text":"","category":"section"},{"location":"develop/","page":"Development","title":"Development","text":"Go to PKGROOT, start a new Julia session and run","category":"page"},{"location":"develop/","page":"Development","title":"Development","text":"julia> using Pkg; Pkg.instantiate()","category":"page"},{"location":"develop/#How-to-build-docs","page":"Development","title":"How to build docs","text":"","category":"section"},{"location":"develop/","page":"Development","title":"Development","text":"Usually, the up-to-state doc is available in here, but there are cases where users need to build the doc themselves.","category":"page"},{"location":"develop/","page":"Development","title":"Development","text":"After instantiating the project, go to PKGROOT, run (without the $ prompt)","category":"page"},{"location":"develop/","page":"Development","title":"Development","text":"$ julia --color=yes docs/make.jl","category":"page"},{"location":"develop/","page":"Development","title":"Development","text":"in your terminal. In a while a folder PKGROOT/docs/build will appear. Open PKGROOT/docs/build/index.html with your favorite browser and have fun!","category":"page"},{"location":"develop/#How-to-run-tests","page":"Development","title":"How to run tests","text":"","category":"section"},{"location":"develop/","page":"Development","title":"Development","text":"After instantiating the project, go to PKGROOT, run (without the $ prompt)","category":"page"},{"location":"develop/","page":"Development","title":"Development","text":"$ julia --color=yes test/runtests.jl","category":"page"},{"location":"develop/","page":"Development","title":"Development","text":"in your terminal.","category":"page"},{"location":"troubleshooting/#Troubleshooting","page":"Troubleshooting","title":"Troubleshooting","text":"","category":"section"},{"location":"troubleshooting/","page":"Troubleshooting","title":"Troubleshooting","text":"This page collects some possible errors you may encounter and trick how to fix them.","category":"page"},{"location":"troubleshooting/","page":"Troubleshooting","title":"Troubleshooting","text":"If you have additional tips, please submit a PR with suggestions.","category":"page"},{"location":"troubleshooting/#Installation-problems","page":"Troubleshooting","title":"Installation problems","text":"","category":"section"},{"location":"troubleshooting/#Cannot-find-the-Julia-executable","page":"Troubleshooting","title":"Cannot find the Julia executable","text":"","category":"section"},{"location":"troubleshooting/","page":"Troubleshooting","title":"Troubleshooting","text":"Make sure you have Julia installed in your environment. Please download the latest stable Julia for your platform. If you are using macOS, the recommended way is to use Homebrew. If you do not want to install Homebrew or you are using other *nix that Julia supports, download the corresponding binaries. And then create a symbolic link /usr/local/bin/julia to the Julia executable. If /usr/local/bin/ is not in your $PATH, modify your .bashrc or .zshrc and export it to your $PATH. Some clusters, like Habanero, Comet already have Julia installed as a module, you may just module load julia to use it. If not, either install by yourself or contact your administrator.","category":"page"},{"location":"troubleshooting/#Loading-settings","page":"Troubleshooting","title":"Loading settings","text":"","category":"section"},{"location":"troubleshooting/","page":"Troubleshooting","title":"Troubleshooting","text":"See \"Loading settings\" for detailed information.","category":"page"},{"location":"troubleshooting/#Loading","page":"Troubleshooting","title":"Loading","text":"","category":"section"},{"location":"troubleshooting/#Why-is-Julia-compiling/loading-modules-so-slow?-What-can-I-do?","page":"Troubleshooting","title":"Why is Julia compiling/loading modules so slow? What can I do?","text":"","category":"section"},{"location":"troubleshooting/","page":"Troubleshooting","title":"Troubleshooting","text":"If you just want Julia do a simple task and only once, you could start Julia REPL with","category":"page"},{"location":"troubleshooting/","page":"Troubleshooting","title":"Troubleshooting","text":"$ julia --compile=min","category":"page"},{"location":"troubleshooting/","page":"Troubleshooting","title":"Troubleshooting","text":"to minimize compilation or","category":"page"},{"location":"troubleshooting/","page":"Troubleshooting","title":"Troubleshooting","text":"$ julia --optimize=0","category":"page"},{"location":"troubleshooting/","page":"Troubleshooting","title":"Troubleshooting","text":"to minimize optimizations, or just use both. Or you could make a system image and run with","category":"page"},{"location":"troubleshooting/","page":"Troubleshooting","title":"Troubleshooting","text":"$ julia --sysimage custom-image.so","category":"page"},{"location":"troubleshooting/","page":"Troubleshooting","title":"Troubleshooting","text":"See Fredrik Ekre's talk for details.","category":"page"},{"location":"install/#Installation","page":"Installation","title":"Installation","text":"","category":"section"},{"location":"install/","page":"Installation","title":"Installation","text":"To install this package, first, you need to install a julia executable from its official website. Version v1.0.0 and above is required. This package may not work on v0.7 and below.","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"If you are using a Mac, and have Homebrew installed, open Terminal.app and type:","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"brew cask install julia","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"Now I am using macOS as a standard platform to explain the following steps:","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"Open Terminal.app, and type julia to start a Julia session.\nRun\njulia> using Pkg; Pkg.update()\n\njulia> pkg\"add Express.jl.git\"\nand wait for its finish.\nRun\njulia> using Express\nand have fun!\nWhile using, please keep this Julia session alive. Restarting may recompile the package and cost some time.","category":"page"},{"location":"install/#Reinstall","page":"Installation","title":"Reinstall","text":"","category":"section"},{"location":"install/","page":"Installation","title":"Installation","text":"In the same Julia session, run\njulia> Pkg.rm(\"Express\"); Pkg.gc()\nPress ctrl+d to quit the current session. Start a new Julia session and repeat the above steps.","category":"page"},{"location":"api/EquationOfStateWorkflow/#Express.EquationOfStateWorkflow-module","page":"EquationOfStateWorkflow module","title":"Express.EquationOfStateWorkflow module","text":"","category":"section"},{"location":"api/EquationOfStateWorkflow/","page":"EquationOfStateWorkflow module","title":"EquationOfStateWorkflow module","text":"Modules = [Express.EquationOfStateWorkflow]","category":"page"},{"location":"configuration/#Configuration-files","page":"Configuration files","title":"Configuration files","text":"","category":"section"},{"location":"configuration/","page":"Configuration files","title":"Configuration files","text":"Express can be run from a configuration file, with some preset rules. The following sections introduce how to write such configuration files. By now, only YAML, JSON, and TOML formats are supported. Please refer their official documentation for their syntax. In the examples below, I will use YAML syntax for configuration files.","category":"page"},{"location":"configuration/#Fitting-equations-of-state","page":"Configuration files","title":"Fitting equations of state","text":"","category":"section"},{"location":"configuration/","page":"Configuration files","title":"Configuration files","text":"Next, specify how many computing cores you are using in total. If running on multiple nodes, write the summation of cores on these nodes. It will be better if np is an integer multiple of the number of pressures.\nThe next key bin is the software being used. If using Quantum ESPRESSO, write qe. The values of qe is the path to the binaries that actually run the calculations. If they are already in the PATH environment variable, write the names of the executables.\ntemplates is a vector of the template files for scf and vc-relax calculations. If it has only one value, all the pressures share the same template. If it has more than one value, the number of files must be equal to the number of pressures. That is, each pressure has its own template.\npressures are the pressures that on the compression curve, they are usually the desired pressures for further calculations. The unit of them is by default GPa. It is usually standard to have at least 6 pressures and at least 1 negative pressure.\ntrial_eos is the starting equation of state for setting volumes for corresponding pressures.\nname is the name of that equation of state. Available options are m (Murnaghan EOS), bm2 - bm4 (Birch–Murnaghan second to fourth order EOSs) and v (Vinet EOS).\nparameters are the parameters of that equation of state. With the first parameter be zero-pressure volume (V_0), the second be zero-pressure bulk modulus (B_0), the third be zero-pressure bulk modulus derivative (B_0). Each value has an associate unit. Allowed units for V_0 are angstrom^3, bohr^3, nm^3, pm^3, etc. Allowed values for B_0 are GPa, Pa, Mbar, kbar, eV/angstrom^3, eV/bohr^3, eV/nm^3, Ry/angstrom^3, Ry/bohr^3, hartree/angstrom^3, etc. B_0 is a dimensionless number so its unit must be 1.\nuse_shell: Whether create shell files to run external software to do computations. Usually this is preferred when express is run non-interactively. If run in interactive mode, you may want to set it to false.","category":"page"},{"location":"configuration/","page":"Configuration files","title":"Configuration files","text":"workflow: eos\nnp: 24\nbin:\n  qe:\n    - pw.x\ntemplates:\n  - examples/Ge/template.in\npressures:\n  unit: GPa\n  values:\n    - -5\n    - -2\n    - 0\n    - 5\n    - 10\n    - 15\n    - 17\n    - 20\ntrial_eos:\n  name: bm3\n  parameters:\n    - - 300.44\n      - bohr^3\n    - - 74.88\n      - GPa\n    - - 4.82\n      - 1\n    - - -612.43149513\n      - Ry\nuse_shell: true","category":"page"},{"location":"configuration/#Phonon-density-of-states-or-phonon-dispersion-relation","page":"Configuration files","title":"Phonon density of states or phonon dispersion relation","text":"","category":"section"},{"location":"configuration/","page":"Configuration files","title":"Configuration files","text":"The first key workflow means which workflow you want to apply. For Phonon density of states, write vdos. For phonon dispersion relation, write phonon dispersion.\nNext, specify how many computing cores you are using in total. If running on multiple nodes, write the summation of cores on these nodes. It will be better if np is an integer multiple of the number of pressures.\nThe next key bin is the software being used. If using Quantum ESPRESSO, write qe. The values of qe is the path to the binaries that actually run the calculations. If they are already in the PATH environment variable, write the names of the executables. They should be written sequentially, with each one corresponds to one calculation.\ntemplates is a vector of the template files for the calculations. If it has only one value, all the pressures share the same template. If it has more than one value, the number of files must be equal to the number of pressures. That is, each pressure has its own template. Note each calculation has a vector of template files.\npressures are the pressures that on the compression curve, they are usually the desired pressures for further calculations. The unit of them is by default GPa. It is usually standard to have at least 6 pressures and at least 1 negative pressure.\nuse_shell: Whether create shell files to run external software to do computations. Usually this is preferred when express is run non-interactively. If run in interactive mode, you may want to set it to false.","category":"page"},{"location":"configuration/","page":"Configuration files","title":"Configuration files","text":"workflow: vdos\nnp: 24\nbin:\n  qe:\n    - pw.x\n    - ph.x\n    - q2r.x\n    - matdyn.x\ntemplates:\n  - - examples/Ge/template.in\n  - - examples/Ge/ph.in\n  - - examples/Ge/q2r.in\n  - - examples/Ge/disp.in\npressures:\n  unit: GPa\n  values:\n    - -5\n    - -2\n    - 0\n    - 5\n    - 10\n    - 15\n    - 17\n    - 20\nuse_shell: true","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = Express","category":"page"},{"location":"#Express","page":"Home","title":"Express","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for Express.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [Express]","category":"page"}]
}
