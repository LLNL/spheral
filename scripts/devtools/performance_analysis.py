"""
Compare performance data for Spheral

Use on LC systems with the following steps:

  1. Run the performance test in Spheral

     $> ./spheral-ats tests/performance.py --logs perftests

  2. Load the virtual environment for Thicket for a bash terminal:

     $> source /usr/gapps/Spheral/venv_timer/bin/activate

     or for a tcsh terminal

     $> source /usr/gapps/Spheral/venv_timer/bin/activate.csh

  3. Run this script and point to the directory created by ATS in step 1

     $> python3 performance_analysis.py --perf-dir perftests
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
def compute_threshold(sf, metric, percent):
    """
    Compute the upper threshold of time using (mu + 2*sigma)(1+percent)
    Parameters:
    sf: Thicket statsframe
    metric: Metric to get upper time limit on
    percent: Additional percent to account for fluctuations
    Returns:
    out: Upper time limit
    """
    mu = gf.dataframe[metric+"_mean"]
    if (metric+"_std" in gf):
        sigma = gf.dataframe[metric+"_std"]
        return (1. + percent)*(mu + 2. * sigma)
    return (1. + percent)*mu

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
    "Remove rows with NANs in a GraphFrame or Thicket"
    if (type(gf) is th.thicket.Thicket):
        query = th.query.Query().match("+", lambda row: row[metric].apply(lambda x: not np.isnan(x)).all())
        return gf.query(query)
    else:
        query = ht.query.Query().match("+", lambda x: not np.isnan(x[metric]))
        return gf.filter(query)

def filter_tests(data1, data2):
    """
    Filters and removes NaNs from input data. Returns two lists of Thickets
    where each entry corresponds to a comparable test. Do not use for scaling
    tests or comparisons between machines, specs, or configurations.
    """
    # Top level metadata comparisons, this should be implicitly handled
    # through benchmark data directory structures
    group1 = ["test_name", "install_config", "cluster"]
    # Lower level metadata comparisons
    group2 = ["numhosts", "jobsize", "total_internal_nodes", "total_steps", "threads_per_rank"]
    filt1 = data1.groupby(group1)
    filt2 = data2.groupby(group1)
    filtered_data1 = []
    filtered_data2 = []
    for test, df1 in filt1.items():
        if test not in filt2:
            print(f"Warning: {test[0]} not found in both data sets, skipping this test")
            continue
        df2 = filt2[test]
        mdata1 = df1.get_unique_metadata()
        mdata2 = df2.get_unique_metadata()
        do_comp = True
        for k in group2:
            val1 = mdata1[k]
            val2 = mdata2[k]
            if (val1 != val2):
                print(f"Warning: {k} does not match for test {test[0]}, {val1} != {val2}.\n"+\
                      "Skipping as tests are not compatible for comparison.")
                do_comp = False
                break
        if (not do_comp):
            continue
        filtered_data1.append(remove_nans(df1))
        filtered_data2.append(remove_nans(df2))
    return filtered_data1, filtered_data2

#---------------------------------------------------------------------------
# Setup argument parser
#---------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--perf-dir", type=str,
                    help="Directory either containing an atsr.py file or a collection of Caliper files")
parser.add_argument("--ref", type=str, default=None,
                    help="Directory of Caliper files to use as reference timings")
parser.add_argument("--display", action="store_true",
                    help="Display a tree for timers that failed")
args = parser.parse_args()
# Create a Thicket of the current performance data
cali_files = []
atsFile = os.path.join(args.perf_dir, "atsr.py")
benchmarks = None
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
    newpath = os.path.join(args.perf_dir, "**/*.cali")
    print(f"Searching {newpath}")
    cali_files = glob.glob(newpath, recursive=True)
if (len(cali_files) == 0):
    raise Exception("No .cali files found")
else:
    print(f"Found {len(cali_files)} Caliper files")

curdata = th.Thicket.from_caliperreader(cali_files)
# Get install config and machine name from current data
install_config = curdata.metadata["install_config"].iloc[0]
machine_name = curdata.metadata["cluster"].iloc[0]
# If no ref or benchmark_dir is provided, look for the benchmark in the
# atsr.py file or a Caliper file
ref_files = args.ref
if (not ref_files):
    # Check in atsr.py
    ref_files = benchmarks
    # Check in a Caliper file
    if (not ref_files):
        try:
            ref_files = curdata.metadata["benchmark_dir"].iloc[0]
        except:
            pass
ref_loc = os.path.join(ref_files,
                       install_config,
                       machine_name,
                       "/latest/*.cali")
cali_ref_files = glob.glob(ref_loc, recursive=True)
if (len(cali_ref_files) == 0):
    # Likely not pointing to benchmark directory so go searching for any Caliper files
    cali_ref_files = glob.glob(os.path.join(ref_files, "**/*.cali"), recursive=True)
if (len(cali_ref_files) > 0):
    refdata = th.Thicket.from_caliperreader(cali_ref_files)
else:
    raise Exception(f"No Caliper files found in {ref_files}")

# In this script, a configuration will refer to the Caliper
# outputs for a specific run of performance.py
# Each configuration is unique to an install config and machine name (cluster)
# and contains at least one Caliper file per test.
# Our current data, provided by the --perf_dir input, should only have
# 1 configuration with 1 Caliper file per test

# Filter both sets of data for tests that are compatible for comparison
cur_config, ref_config = filter_tests(curdata, refdata)
tests_failed = []
# Iterate over each test in the current config
for cprof, cref in zip(cur_config, ref_config):
    mult_refs = True
    # Get statistical data for reference config
    if (len(cref.profile) == 1):
        print(f"Warning: Only 1 reference run found for {test}")
        mult_refs = False
    metric0 = "Avg time/rank"
    metric1 = "Avg time/rank (exc)"
    metrics = [metric0, metric1]
    # Get statistical values
    th.stats.mean(cref, metrics)
    if mult_refs:
        th.stats.std(cref, metrics)
    # Compute the max allowable time for the main region
    cprof.statsframe.dataframe["upper_limit"] = compute_upper_limit(cref.statsframe, metric0, percent)
    cprof.move_metrics_to_statsframe([metric0])
    cprof.statsframe.dataframe["diff"] = cprof.statsframe.dataframe["upper_limit"] - cprof.statsframe.dataframe[metric0]
    # See if the time of "main" exceeds the upper limit
    cur_diff = get_times(cprof.statsframe, "main", "diff")
    if (cur_diff < 0.):
        test_failed.append(cprof)
        print(f"'main' in test {test} has exceeded the upper limit")
        if args.display:
            display(cprof.statsframe.tree(metric1, metric1+"_mean"))
if (len(test_failed) > 0):
    raise Exception(f"{len(test_failed)} have failed")
