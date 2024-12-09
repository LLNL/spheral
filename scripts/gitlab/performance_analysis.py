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

     $> python3 performance_analysis.py --files perftests
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

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# Location of benchmark data
benchmark_dir = "/usr/gapps/Spheral/benchmarks"

class PerfProcess:
    def __init__(self):
        #---------------------------------------------------------------------------
        # Setup argument parser
        #---------------------------------------------------------------------------
        parser = argparse.ArgumentParser()
        parser.add_argument("--files", type=str, required=True,
                            help="Directory either containing an atsr.py file or a collection of Caliper files")
        parser.add_argument("--ref-files", type=str, default=benchmark_dir,
                            help="Directory of Caliper files to use as reference timings")
        args = parser.parse_args()
        cali_files = self.get_cali_files(args.files)
        if (len(cali_files) == 0):
            raise Exception("No .cali files found")
        else:
            print(f"Found {len(cali_files)} Caliper files")
        self.newdata = th.Thicket.from_caliperreader(cali_files)
        cali_ref_files = self.get_cali_files(args.ref_files)
        if (len(cali_ref_files) > 0):
            self.refdata = th.Thicket.from_caliperreader(cali_ref_files)
        else:
            self.refdata = None

    def get_cali_files(self, lfiles):
        "Retrieve Caliper files either using an atsr.py file or just glob"
        cali_files = []
        atsFile = os.path.join(lfiles, "atsr.py")
        if (os.path.exists(atsFile)):
            # Run atsr.py and put values into globals
            exec(compile(open(atsFile).read(), atsFile, 'exec'), globals())
            state = globals()["state"]
            tests = [t for t in state["testlist"] if t['status'] == PASSED]
            for test in tests:
                # Retrieve the Caliper file from run
                run_dir = test["directory"]
                cali_file = test["options"]["caliper_filename"]
                cfile = os.path.join(run_dir, cali_file)
                cali_files.append(cfile)
                print(cfile)
        else:
            newpath = os.path.join(lfiles, "**/*.cali")
            print(f"Searching {newpath}")
            cali_files = glob.glob(newpath, recursive=True)
        return cali_files

    def filter_tests(self, thicket):
        """
        Return a dictionary of Thickets where each entry is for a different test
        """
        return thicket.groupby(["test_name", "install_config", "numhosts"])

    def remove_nans(self, gf, metric):
        "Remove rows with NANs in a GraphFrame or Thicket"
        if (type(gf) == th.thicket.Thicket):
            query = th.query.Query().match("+", lambda row: row[metric].apply(lambda x: x >= 0).all())
            return gf.query(query)
        elif:
            query = ht.query.Query().match("+", lambda x: x[metric] >= 0.])
            return gf.filter(query)

    def compare_times(self, metric="Avg time/rank (exc)", region="main"):
        "Compare times between the new data and the reference data"
        if not self.refdata:
            raise Exception("Cannot compare times, no reference data found")
        threfs = self.filter_tests(self.refdata)
        thtests = self.filter_tests(self.newdata)
        # Loop over the tests/install/num hosts configs
        for test, cprof in thtests.items():
            # Grab reference data or skip if nothing to compare to
            if test in threfs:
                cref = threfs[test]
            else:
                continue
            th.stats.mean(cref, [metric])
            th.stats.std(cref, [metric])
            cref.statsframe = self.remove_nans(cref.statsframe, metric+"_mean")
            cref.statsframe = self.remove_nans(cref.statsframe, metric+"_std")
            statdf = cref.statsframe
            # Ensure we only have 1 profile (.cali file) per test
            if (len(cprof.profile) > 1):
                raise Exception(f"Multiple profiles found for {test}")
            filtth = self.remove_nans(cprof, metric)
            filtth.move_metrics_to_statsframe(metric)
            # Find the relative difference between the single test and the mean of the reference tests
            filtth.statsframe.dataframe["diff"] = (filtth.statsframe.dataframe[metric] -
                                                   statdf.dataframe[metric+"_mean"])/statdf.dataframe[metric+"_mean"]

sphp = PerfProcess()
sphp.compare_times()
