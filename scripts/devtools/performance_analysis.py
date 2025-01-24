"""
Compare performance data for Spheral

Use on LC systems with the following steps:

  1. Run the performance test in Spheral

     $> ./spheral-ats --numNodes 2 --logs test_dir_name tests/performance.py

  2. Load the virtual environment for Thicket for a bash terminal:

     $> source /usr/gapps/Spheral/venv_timer/bin/activate

     or for a tcsh terminal

     $> source /usr/gapps/Spheral/venv_timer/bin/activate.csh

  3. Run this script and point to the directory created by ATS in step 1

     $> python3 performance_analysis.py --perf-dir test_dir_name
"""

import os, sys, shutil, glob
import argparse

try:
    import thicket as th
    import hatchet as ht
except:
    print("Thicket not found. Be sure to load virtual environment first")
    raise Exception

from IPython.display import display
from IPython.display import HTML

import numpy as np

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

percent = 0.08
def compute_threshold(sf, metric):
    """
    Compute the threshold of time as thresh = 0.08*mu + 2*sigma
    Will be used in the comparison abs(t - mu) > thresh to determine
    notable differences in times
    Parameters:
    sf: Thicket statsframe (Hatchet Graphframe) for reference data
    metric: Metric for threshold
    Returns:
    thresh: Time limit threshold Graphframe variable
    """
    mu = sf.dataframe[metric+"_mean"]
    if (metric+"_std" in sf.dataframe):
        sigma = sf.dataframe[metric+"_std"]
        return percent*mu + 2.*sigma
    return percent*mu

def get_times(gf, region, metric = "Avg time/rank"):
    """
    Return the times for a given region and metric.
    If multiple regions have the same name, the sum of those times are returned.
    Returns an array of times for each profile in the Thicket
    """
    if (type(gf) is th.thicket.Thicket):
        cregion = gf.get_node(region)
        profs = gf.profile
        if len(profs) > 1:
            times = []
            for i in profs:
                vth = gf.filter_profile([i])
                times.append(sum(vth.dataframe.loc[cregion][metric].values))
        else:
            times = [sum(gf.dataframe.loc[cregion][metric].values)]
    else:
        names = gf.dataframe["name"]
        indx = names[names == region].index
        times = [sum(gf.dataframe.loc[indx][metric].values)]
    return times

def remove_nans(gf, metric="Avg time/rank"):
    "Remove rows with NANs in a GraphFrame or Thicket or list/dict of those types"
    if (type(gf) is dict or type(gf) is th.groupby.GroupBy):
        newdict = {}
        for key, val in gf.items():
            newval = remove_nans(val)
            newdict.update({key: newval})
        return newdict
    elif (type(gf) is list):
        newlist = []
        for val in gf:
            newval = remove_nans(val)
            newlist.append(newval)
        return newlist
    elif (type(gf) is th.thicket.Thicket):
        query = th.query.Query().match("+", lambda row: row[metric].apply(lambda x: not np.isnan(x)).all())
        return gf.query(query)
    elif (type(gf) is ht.GraphFrame):
        query = ht.query.Query().match("+", lambda x: not np.isnan(x[metric]))
        return gf.filter(query)
    else:
        raise TypeError(f"Unrecognized type in remove_nans {type(gf)}")

def group_tests(data):
    """
    Groups input data based on tests and removes NaNs.
    Parameters:
    data: Thicket to filter
    Returns:
    filt: GroupBy of Thickets based on tests
    """
    test_group = ["test_name", "total_internal_nodes", "total_steps"]
    filt = data.groupby(test_group)
    return remove_nans(filt)

def compare_metadata(cdata, rdata, tests):
    cmdata = cdata.get_unique_metadata()
    rmdata = rdata.get_unique_metadata()
    failed_configs = []
    for t in tests:
        cval = cmdata[t]
        rval = rmdata[t]
        if (cval != rval):
            failed_configs.append(f"{t}: {cval} != {rval}")
    return failed_configs

def compare_tests(cdata, rdata):
    "Determine why a test configuration differs"
    tests = ["total_internal_nodes", "total_steps"]
    return compare_metadata(cdata, rdata, tests)

def compare_config(cdata, rdata):
    """
    Check if tests between two sets of data are compatible for direct comparison.
    Do not use if comparing across hardware or install configurations.
    Assumes Thicket data has comparable test_name, total_internal_nodes, and total_steps
    """
    hardware_tests = ["install_config", "cluster", "jobsize", "threads_per_rank"]
    return compare_metadata(cdata, rdata, hardware_tests)

def filter_tests(data, test_name):
    return data.filter_metadata(lambda x: x["test_name"] == test_name)

def get_caliper_files(file_path):
    atsFile = os.path.join(file_path, "atsr.py")
    cali_files = []
    benchmarks = None
    # If perf-dir is an ATS output directory, find the Caliper files from atsr.py
    if (os.path.exists(atsFile)):
        # Run atsr.py and put values into globals
        exec(compile(open(atsFile).read(), atsFile, 'exec'), globals())
        state = globals()["state"]
        tests = [t for t in state["testlist"] if t['status'] == PASSED]
        for test in tests:
            # Retrieve the Caliper file from run
            run_dir = test["directory"]
            cali_file = test["options"]["caliper_filename"]
            # Check if benchmark_dir is in ats options
            if (not benchmarks and "benchmark_dir" in test["options"]):
                benchmarks = test["options"]["benchmark_dir"]
            cfile = os.path.join(run_dir, cali_file)
            cali_files.append(cfile)
    else:
        newpath = os.path.join(file_path, "**/*.cali")
        print(f"Searching {newpath}")
        cali_files = glob.glob(newpath, recursive=True)
    return cali_files, benchmarks

#---------------------------------------------------------------------------
# Setup argument parser
#---------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--perf-dir", type=str,
                    help="Directory either containing an atsr.py file or a collection of Caliper files")
parser.add_argument("--ref", type=str, default=None,
                    help="Directory of Caliper files to use as reference timings.\n"+\
                    "Uses the shared benchmark performance data from CI by default.")
parser.add_argument("--diff-configs", action="store_true",
                    help="Set to true to allow comparisons across hardware/compiler configuration")
parser.add_argument("--display", action="store_true",
                    help="Display a tree for timers that failed")
args = parser.parse_args()

# Create a Thicket of the current performance data
#-------------------------------------------------

cali_files, benchmarks = get_caliper_files(args.perf_dir)
if (len(cali_files) == 0):
    raise Exception(f"No .cali files found in {args.perf_dir}")
curdata = th.Thicket.from_caliperreader(cali_files)

# Create a Thicket of the reference data
#---------------------------------------

# If reference data is specified, grab the data directly
if (args.ref):
    cali_ref_files, unused_benchmark = get_caliper_files(args.ref)
else:
    # If no ref or benchmark_dir is provided, look for the benchmark in the
    # atsr.py file or a Caliper file
    # Check for reference from atsr.py
    ref_files = benchmarks
    # Check in a Caliper file
    if (not ref_files):
        try:
            ref_files = curdata.metadata["benchmark_dir"].iloc[0]
        except:
            raise Exception("No reference or benchmark data specified")
    # If we using benchmark reference data, only grab the current install/machine
    # Get install config and machine name from current data
    install_config = curdata.metadata["install_config"].iloc[0]
    machine_name = curdata.metadata["cluster"].iloc[0]
    ref_loc = os.path.join(ref_files,
                           install_config,
                           machine_name,
                           "latest/*.cali")
    cali_ref_files = glob.glob(ref_loc, recursive=True)
if (len(cali_ref_files) == 0):
    raise Exception(f"No Caliper files found in {cali_ref_files}")
refdata = th.Thicket.from_caliperreader(cali_ref_files)

# Group, filter, and compare performance data
#--------------------------------------------

# Filter both sets of data set by the tests
cur_test_data = group_tests(curdata)
ref_test_data = group_tests(refdata)

test_status = {}
failed_tests = {}
for test_key, ctest in cur_test_data.items():
    test_name = test_key[0]
    test_sph_nodes = test_key[1]
    test_steps = test_key[2]
    if (test_key not in ref_test_data):
        # This means the test configurations differ (number of time steps etc)
        if (test_name not in refdata.get_unique_metadata()["test_name"]):
            skip_msg = f"{test_name} not found in reference data"
            test_status.update({test_name: ("SKIPPED-TEST", [skip_msg])})
        else:
            rtest = filter_tests(refdata, test_name)
            ftest_configs = compare_tests(ctest, rtest)
            test_status.update({test_name: ("SKIPPED-TEST", ftest_configs)})
        continue
    rtest = ref_test_data[test_key]
    if (not args.diff_configs):
        fh_configs = compare_config(ctest, rtest)
        if (fh_configs):
            # This means the hardware/compiler configurations differs
            test_status.update({test_name: ("SKIPPED-CONF", fh_configs)})
            continue
    mult_refs = True
    # Get statistical data for reference config
    if (len(rtest.profile) == 1):
        print(f"Warning: Only 1 reference run found for {test_name}")
        mult_refs = False
    metric0 = "Avg time/rank"
    metric1 = "Avg time/rank (exc)"
    metrics = [metric0, metric1]
    # Get statistical values
    th.stats.mean(rtest, metrics)
    if mult_refs:
        th.stats.std(rtest, metrics)
    # Compute the max allowable time for the main region
    ctest.statsframe.dataframe["thresh"] = compute_threshold(rtest.statsframe, metric0)
    ctest.move_metrics_to_statsframe([metric0])
    ref_main = get_times(rtest.statsframe, "main", metric0+"_mean")[0]
    cur_main = get_times(ctest.statsframe, "main", metric0)[0]
    ref_thresh = get_times(ctest.statsframe, "main", "thresh")[0]
    main_diff = cur_main - ref_main
    if (main_diff > ref_thresh):
        cur_status = "FAILED"
        if args.display:
            display(ctest.statsframe.tree(metric1, metric1+"_mean"))
    elif (main_diff < -ref_thresh):
        cur_status = "PASSED"
        if args.display:
            display(ctest.statsframe.tree(metric1, metric1+"_mean"))
    else:
        cur_status = "PASSED"
    test_status.update({test_name: (cur_status, cur_main, ref_main, ref_thresh)})
num_failed = 0
for test_name, val in test_status.items():
    if ("SKIPPED" in val[0]):
        print(f"{test_name}: SKIPPED. Differences found for:")
        for i in val[1]:
            print(i)
        if (val[0] == "SKIPPED-CONF"):
            print("Rerun with --diff-configs to allow comparisons across hardware/compilers")
    else:
        ctime = val[1]
        rtime = val[2]
        thresh = val[3]
        if ("FAILED" in val[0]):
            num_failed += 1
            print(f"{test_name}: FAILED. Time of 'main' current/ref: {ctime/rtime*100.:0.3f}%")
        else:
            print(f"{test_name}: PASSED. Time of 'main' current/ref: {ctime/rtime*100.:0.3f}%")
if (num_failed > 0):
    raise Exception(f"{num_failed} have failed")
