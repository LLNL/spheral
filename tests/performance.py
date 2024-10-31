#!/user/bin/env python3

# This file runs and compares performance tests through the ats system.
# Run using: ./spheral-ats tests/performance.py

import sys, shutil, os, time
import numpy as np
cur_dir = os.path.dirname(__file__)
spheral_path = os.path.join(cur_dir, "../lib/python3.9/site-packages/Spheral")
sys.path.append(spheral_path)
import SpheralConfigs

# Current system architecture from Spack
spheral_sys_arch = SpheralConfigs.sys_arch()
# Current install configuration from Spack
spheral_install_config = SpheralConfigs.config()

glue(keep=True)

def add_timer_cmds(cali_name, test_name):
    return f"--caliperFilename {cali_name} --adiakData 'test_name: {test_name}, install_config: {spheral_install_config}'"

if ("power" in spheral_sys_arch):
    num_nodes = 1
    num_cores = 40
elif ("broadwell" in spheral_sys_arch):
    num_nodes = 2
    num_cores = 36

# NOH tests
test_dir = os.path.join(SpheralConfigs.test_install_path(), "functional/Hydro/Noh")

# Select which timing regions to post-process
regions = ["CheapRK2",
           "CheapRK2PreInit",
           "ConnectivityMap_computeConnectivity",
           "ConnectivityMap_patch",
           "CheapRK2EvalDerivs",
           "CheapRK2EndStep"]
# Select which timers to use to post-process the regions above
timers = ["sum#inclusive#sum#time.duration"] # Means the sum of the time from all ranks

# General input for all Noh tests
gen_noh_inps = "--crksph False --cfl 0.25 --Cl 1.0 --Cq 1.0 --xfilter 0.0 "+\
    "--nPerh 2.01 --graphics False --clearDirectories False --doCompare False "+\
    "--dataDir None --vizTime None --vizCycle None"

group(name="NOH 2D tests")
test_file = "Noh-cylindrical-2d.py"
nRadial = 100
test_path = os.path.join(test_dir, test_file)
test_name_base = "NC2D"

# Test with varying number of ranks
ranks = [1, 2, 4]
for i, n in enumerate(ranks):
    test_name = f"{test_name_base}_{i}"
    caliper_filename = f"{test_name}_{int(time.time())}.cali"
    timer_cmds = add_timer_cmds(caliper_filename, test_name)
    inps = f"{gen_noh_inps} --nRadial {nRadial} --steps 10 {timer_cmds}"
    ncores = int(num_nodes*num_cores/n)
    t = test(script=test_path, clas=inps, label=f"{test_name}",
             np=ncores,
             caliper_filename=caliper_filename,
             regions=regions,
             timers=timers,
             install_config=spheral_install_config)

endgroup()

group(name="NOH 3D tests")
test_file = "Noh-spherical-3d.py"
test_path = os.path.join(test_dir, test_file)
test_name_base = "NS3D"

# Test with varying number of SPH nodes per rank
npcore = [100, 200, 300]
for i, n in enumerate(npcore):
    test_name = f"{test_name_base}_{i}"
    caliper_filename = f"{test_name}_{int(time.time())}.cali"
    total_sph_nodes = n*num_cores
    npd = int(np.cbrt(total_sph_nodes))
    node_inps = f"--nx {npd} --ny {npd} --nz {npd}"
    timer_cmds = add_timer_cmds(caliper_filename, test_name)
    inps = f"{gen_noh_inps} {node_inps} --steps 3 {timer_cmds}"
    # WIP: Path to benchmark timing data
    ncores = int(num_cores)
    t = test(script=test_path, clas=inps, label=f"{test_name}",
             np=ncores,
             independent=False,
             caliper_filename=caliper_filename,
             regions=regions,
             timers=timers,
             install_config=spheral_install_config)
# Add a wait to ensure all timer files are done
wait()
