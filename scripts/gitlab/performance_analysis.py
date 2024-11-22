#!/user/bin/env python3

import os, sys, shutil, glob
import argparse

import SpheralConfigs

# Location of benchmark data
benchmark_dir = "/usr/gapps/Spheral/benchmarks"

caliper_loc = SpheralConfigs.caliper_module_path()
sys.path.append(caliper_loc)
import caliperreader as cr

def main():
    #---------------------------------------------------------------------------
    # Setup argument parser
    #---------------------------------------------------------------------------
    parser = argparse.ArgumentParser()
    parser.add_argument("--atsOutput", type=str, required=True,
                        help="Path to atsr.py file produced from running performance.py")
    args = parser.parse_args()

    atsFile = args.atsOutput
    if (os.path.isdir(args.atsOutput)):
        atsFile = os.path.join(args.atsOutput, "atsr.py")
    if (not os.path.exists(atsFile)):
        raise Exception("ATS file not found")
    # Run atsr.py and put values into globals
    exec(compile(open(atsFile).read(), atsFile, 'exec'), globals())
    state = globals()["state"]
    tests = state["testlist"]
    for test in tests:
        # Retrieve the Caliper file from run
        run_dir = test["directory"]
        options = test["options"]
        cali_file = options["caliper_filename"]
        cfile = os.path.join(run_dir, cali_file)
        install_config = options["install_config"]

        # Grab list of regions and timers
        ref_regions = options["regions"]
        ref_timers = options["timers"]
        # Read Caliper file
        r = cr.CaliperReader()
        r.read(cfile)
        # Get adiak metadata
        gls = r.globals
        test_name = gls["test_name"]

        # Extract relevant regions and timers
        times = {}
        for rec in records:
            if ("region" in rec):
                fname = rec["region"]
                if (type(fname) is list):
                    fname = fname[-1]
                if (fname in ref_regions):
                    if (fname in times):
                        for t in ref_timers:
                            times[fname][t] += float(rec[t])
                    else:
                        new_dict = {}
                        for t in ref_timers:
                            new_dict.update({t: float(rec[t])})
                        times.update({fname: new_dict})
        # Get historical timing data
        cali_ref_dir = os.path.join(benchmark_dir, install_config, test_name)
        if (not os.path.exists(cali_ref_dir)):
            os.makedirs(cali_ref_dir)
        shutils.copyfile(cfile, os.path.join(cali_ref_dir, cali_file))

if __name__=="__main__":
    main()
