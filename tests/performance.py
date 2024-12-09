#!/user/bin/env python3

# This file runs and compares performance tests through the ats system.
# Run using: ./spheral-ats tests/performance.py

import sys, shutil, os, time, stat
import numpy as np
spheral_path = "../lib/python3.9/site-packages/Spheral"
sys.path.append(spheral_path)
import SpheralConfigs
sys.path.append("../scripts")
from spheralutils import num_3d_cyl_nodes
from spheralutils import num_2d_cyl_nodes
from ats import configuration

# Get options from ats
opts = getOptions()

# For CI runs, automatically copy caliper files to benchmark location
output_loc = None
if "cirun" in opts and opts["cirun"]:
    output_loc = SpheralConfigs.benchmark_data()

# Current system architecture from Spack
spheral_sys_arch = SpheralConfigs.sys_arch()
# Current install configuration from Spack
spheral_install_config = SpheralConfigs.config()

# Retrieve the host name and remove any numbers
temp_uname = os.uname()
hostname = "".join([i for i in temp_uname[1] if not i.isdigit()])
mac_procs = {"rzhound": 112, "rzwhippet": 112, "ruby": 112,
             "rzadams": 96, "rzvernal": 64, "tioga": 64,
             "rzansel": 40, "lassen": 40, "rzgenie": 36}
# Find out how many nodes our allocation has grabbed
num_nodes = max(1, configuration.machine.numNodes)
num_cores = 0
try:
    num_cores = mac_procs[hostname] * num_nodes
except:
    print("Machine name not recognized")
    raise Exception

# General number of SPH nodes per core
n_per_core_3d = 8000
n_per_core_2d = 1000
# Do not test SVPH
# Test all others with and without solid if possible
# 5k-10k nodes per core for 3d, 1k nodes per core for 2d

def gather_files(manager):
    '''
    Function to gather Caliper file when ATS is finished running
    '''
    filtered = [test for test in manager.testlist if test.status is PASSED]
    # Set read/write permissions for owner and read/write permissions for group
    perms = stat.S_IRUSR | stat.S_IWUSR | stat.S_IWGRP | stat.S_IRGRP
    for test in filtered:
        run_dir = test.directory
        cali_filename = test.options["caliper_filename"]
        cfile = os.path.join(run_dir, cali_filename)
        test_name = test.options["label"]
        outdir = os.path.join(output_loc, spheral_install_config)
        if (not os.path.exists(outdir)):
            log(f"Creating {outdir}")
            os.makedirs(outdir)
        outfile = os.path.join(outdir, cali_filename)
        log(f"Copying {cali_filename} to {outdir}")
        shutil.copy(cfile, outfile)
        os.chmod(outfile, perms)

def spheral_setup_test(test_path, test_name, inps, ncores, threads=1, **kwargs):
    '''
    General method for creating an individual performance test
    '''
    cali_name = f"{test_name}_{int(time.time())}.cali"
    timer_cmds = f"--caliperFilename {cali_name} --adiakData 'test_name: {test_name}'"
    finps = f"{inps} {timer_cmds}"
    t = test(script=test_path, clas=finps,
             label=test_name,
             np=ncores,
             nt=threads,
             caliper_filename=cali_name,
             **kwargs)
    return t

if (output_loc):
    onExit(gather_files)
glue(keep=True)

# Compute number of SPH nodes
Ntotal = num_cores*n_per_core_3d

#---------------------------------------------------------------------------
# Taylor impact test
#---------------------------------------------------------------------------
test_dir = os.path.join(SpheralConfigs.test_install_path(), "functional/Strength/TaylorImpact")

group(name="Taylor impact tests")
test_file = "TaylorImpact.py"
test_path = os.path.join(test_dir, test_file)
test_name = "3DTAYLOR"

# Estimate nr and nz for cylindrical node distribution to have Ntotal nodes
rlen = 0.945
zlen = 7.5
steps = 10
nr, nz = num_3d_cyl_nodes(0., rlen, 0., zlen, 0., 2.*np.pi, 10, 80, Ntotal)
gen_inps = f"--geometry 3d --steps {steps} --compatibleEnergy False "+\
    "--densityUpdate SumVoronoiCellDensity --clearDirectories False --baseDir None "+\
    "--vizTime None --vizCycle None --siloSnapShotFile None "+\
    f"--rlength {rlen} --zlength {zlen} --nr {nr} --nz {nz}"

# Test variations
test_inp = {"CRK": "--crksph True", "FSI": "--fsisph True"}
for tname, tinp in test_inp.items():
    inps = f"{gen_inps} {tinp}"
    t = spheral_setup_test(test_path, test_name+tname, inps, num_cores)

#---------------------------------------------------------------------------
# 3D convection test
#---------------------------------------------------------------------------
test_dir = os.path.join(SpheralConfigs.test_install_path(), "unit/Boundary")

group(name="Convection tests")
test_file = "testPeriodicBoundary-3d.py"
test_path = os.path.join(test_dir, test_file)
test_name = "3DCONV"

# Test with 1 and 2 threads
npd = int(np.cbrt(Ntotal))
steps = 10
inps = f"--nx {npd} --ny {npd} --nz {npd} --steps {steps}"
thrs = [1, 2]
for thr in thrs:
    ncores = int(num_cores / thr)
    t = spheral_setup_test(test_path, test_name+f"THR{thr}", inps, ncores, thr)

#---------------------------------------------------------------------------
# NOH tests
#---------------------------------------------------------------------------
fluid_variations = {"SPH": "--crksph False --solid True",
                    "FSISPH": "--fsisph True --solid True",
                    "CRKSPH": "--crksph True --solid True",
                    "PSPH": "--psph True",
                    "GSPH": "--gsph True",
                    "MFM": "--mfm True",
                    "MFV": "--mfv True"}
test_dir = os.path.join(SpheralConfigs.test_install_path(), "functional/Hydro/Noh")

# General input for all Noh tests
gen_noh_inps = "--cfl 0.25 --Cl 1.0 --Cq 1.0 --xfilter 0.0 "+\
    "--nPerh 2.01 --graphics False --clearDirectories False --doCompare False "+\
    "--dataDir None --vizTime None --vizCycle None"

#++++++++++++++++++++
# Noh cylindrical 2d
#++++++++++++++++++++
group(name="2D NOH tests")
test_file = "Noh-cylindrical-2d.py"
test_path = os.path.join(test_dir, test_file)
test_name = "NC2D"

steps = 100
rmin = 0.
rmax = 1.
thetaFactor = 0.5
# Only use half the number of cores
ncores = int(num_cores/2)
Ntotal2d = n_per_core_2d*ncores
nRadial = num_2d_cyl_nodes(rmin, rmax, thetaFactor*np.pi, 100, Ntotal2d)
gen_inps = f"{gen_noh_inps} --nRadial {nRadial} --steps {steps}"

# Test with different hydro options
for tname, tinp in fluid_variations.items():
    inps = gen_inps + " " + tinp
    t = spheral_setup_test(test_path, test_name+tname, inps, ncores, independent=True)

#++++++++++++++++++++
# Noh spherical 3d
#++++++++++++++++++++
group(name="3D NOH tests")
test_file = "Noh-spherical-3d.py"
test_path = os.path.join(test_dir, test_file)
test_name = "NS3D"

steps = 10
gen_inps = f"{gen_noh_inps} --nx {npd} --ny {npd} --nz {npd} --steps {steps}"
# Test with different hydro options
for tname, tinp in fluid_variations.items():
    inps = gen_inps + " " + tinp
    t = spheral_setup_test(test_path, test_name+tname, inps, num_cores)
# Add a wait to ensure all timer files are done
wait()
