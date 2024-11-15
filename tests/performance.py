#!/user/bin/env python3

# This file runs and compares performance tests through the ats system.
# Run using: ./spheral-ats tests/performance.py

import sys, shutil, os, time
import numpy as np
spheral_path = "../lib/python3.9/site-packages/Spheral"
sys.path.append(spheral_path)
import SpheralConfigs

# If desired, set a location to consolidate Caliper files, tthis is useful
# when running scaling tests
# This automatically creates directories based on the install configuration
# and test names inside output_loc
# WARNING: Be sure to remove older performance data in
# output location beforehand
#output_loc = "/home/user/scaling/test/files"
output_loc = None

# Current system architecture from Spack
spheral_sys_arch = SpheralConfigs.sys_arch()
# Current install configuration from Spack
spheral_install_config = SpheralConfigs.config()

# Consolidate Caliper files after run
def gather_files(manager):
    filtered = [test for test in manager.testlist if test.status is PASSED]
    for test in filtered:
        run_dir = test.directory
        cali_filename = test.options["caliper_filename"]
        cfile = os.path.join(run_dir, cali_filename)
        test_name = test.options["label"]
        outdir = os.path.join(output_loc, spheral_install_config, test_name)
        if (not os.path.exists(outdir)):
            log(f"Creating {outdir}")
            os.mkdir(outdir)
        outfile = os.path.join(outdir, cali_filename)
        log(f"Copying {cali_filename} to {outdir}")
        shutil.copy(cfile, outfile)

if (output_loc):
    onExit(gather_files)
glue(keep=True)

def add_timer_cmds(cali_name, test_name):
    return f"--caliperFilename {cali_name} --adiakData 'test_name: {test_name}, install_config: {spheral_install_config}'"

if ("power" in spheral_sys_arch):
    num_nodes = 1
    num_cores = 40
elif ("broadwell" in spheral_sys_arch):
    num_nodes = 2
    num_cores = 36

# Select which timing regions to compare (for CI)
regions = ["CheapRK2",
           "CheapRK2PreInit",
           "ConnectivityMap_computeConnectivity",
           "ConnectivityMap_patch",
           "CheapRK2EvalDerivs",
           "CheapRK2EndStep"]
# Select which timers to compare (for CI)
timers = ["sum#inclusive#sum#time.duration"] # Means the sum of the time from all ranks

# 3D convection test
test_dir = os.path.join(SpheralConfigs.test_install_path(), "unit/Boundary")

group(name="3D Convection test")
test_file = "testPeriodicBoundary-3d.py"
test_path = os.path.join(test_dir, test_file)
test_name = "3DCONV"

# Test with varying number of ranks
ranks = [1, 2, 4]
# We want 20 points per unit length
ref_len = 1.
sph_point_rho = 20. / ref_len
sph_per_core = 300
for i, n in enumerate(ranks):
    caliper_filename = f"{test_name}_{i}_{int(time.time())}.cali"
    timer_cmds = add_timer_cmds(caliper_filename, test_name)
    ncores = int(num_nodes*num_cores/n)
    total_sph_nodes = sph_per_core * ncores
    npd = int(np.cbrt(total_sph_nodes))
    new_len = npd * ref_len / sph_point_rho
    inps = f"--nx {npd} --ny {npd} --nz {npd} --x1 {new_len} --y1 {new_len} --z1 {new_len} --steps 100 {timer_cmds}"
    t = test(script=test_path, clas=inps,
             label=test_name,
             np=ncores,
             caliper_filename=caliper_filename,
             regions=regions,
             timers=timers,
             install_config=spheral_install_config)

endgroup()

# NOH tests
test_dir = os.path.join(SpheralConfigs.test_install_path(), "functional/Hydro/Noh")

# General input for all Noh tests
gen_noh_inps = "--crksph False --cfl 0.25 --Cl 1.0 --Cq 1.0 --xfilter 0.0 "+\
    "--nPerh 2.01 --graphics False --clearDirectories False --doCompare False "+\
    "--dataDir None --vizTime None --vizCycle None"

group(name="NOH 2D tests")
test_file = "Noh-cylindrical-2d.py"
nRadial = 100
test_path = os.path.join(test_dir, test_file)
test_name = "NC2D"

# Test with varying number of ranks
ranks = [1, 2, 4]
for i, n in enumerate(ranks):
    caliper_filename = f"{test_name}_{i}_{int(time.time())}.cali"
    timer_cmds = add_timer_cmds(caliper_filename, test_name)
    inps = f"{gen_noh_inps} --nRadial {nRadial} --steps 10 {timer_cmds}"
    ncores = int(num_nodes*num_cores/n)
    t = test(script=test_path, clas=inps,
             label=test_name,
             np=ncores,
             caliper_filename=caliper_filename,
             regions=regions,
             timers=timers,
             install_config=spheral_install_config)

endgroup()

group(name="NOH 3D tests")
test_file = "Noh-spherical-3d.py"
test_path = os.path.join(test_dir, test_file)
test_name = "NS3D"

# Test with varying number of SPH nodes per rank
npcore = [100, 200, 300]
for i, n in enumerate(npcore):
    caliper_filename = f"{test_name}_{i}_{int(time.time())}.cali"
    ncores = int(num_nodes*num_cores)
    total_sph_nodes = n*ncores
    npd = int(np.cbrt(total_sph_nodes))
    node_inps = f"--nx {npd} --ny {npd} --nz {npd}"
    timer_cmds = add_timer_cmds(caliper_filename, test_name)
    inps = f"{gen_noh_inps} {node_inps} --steps 10 {timer_cmds}"
    # WIP: Path to benchmark timing data
    t = test(script=test_path, clas=inps,
             label=test_name,
             np=ncores,
             caliper_filename=caliper_filename,
             regions=regions,
             timers=timers,
             install_config=spheral_install_config)
# Add a wait to ensure all timer files are done
wait()
