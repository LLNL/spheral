#-------------------------------------------------------------------------------
# Create a standard and hopefully convenient command line parser for Spheral
# scripts.
#-------------------------------------------------------------------------------
import argparse, mpi

from SpheralCompiledPackages import *

from SpheralTestUtilities import globalFrame
from SpheralUtilities import TimerMgr

def commandLine(**options):

    # Build a command line parser with the keyword arguments passed to us.
    parser = argparse.ArgumentParser()
    for key in options:
        parser.add_argument("--" + key,
                            dest = key,
                            default = options[key])

    # Add the universal options supported by all Spheral++ scripts.
    parser.add_argument("-v", "--verbose",
                        action = "store_true",
                        dest = "verbose",
                        default = False,
                        help = "Verbose output -- print all options that were set.")
    # This logic checks if the user already set a Caliper argument and default value
    # and prevents adding the argument if it already exists
    arg_list = [action.dest for action in parser._actions]
    cali_args = ["Config", "Filename", "ConfigJSON"]
    for ca in cali_args:
        if (ca not in arg_list):
            parser.add_argument(f"--caliper{ca}", default="", type=str)
    # Evaluate the command line.
    args = parser.parse_args()
    arg_dict = vars(args)

    if (not TimerMgr.timers_usable()):
        if (args.caliperConfig or args.caliperFilename or args.caliperConfigJSON):
            print("WARNING: Caliper command line inputs provided for "+\
                  "non-timer install. Reconfigure the install with "+\
                  "-DENABLE_TIMER=ON to be able to use Caliper timers.")

    # Verbose output?
    if args.verbose:
        print("All parameters set:")
        for key, val in arg_dict.items():
            if key in options:
                if val != options[key]:
                    print("  *  ", key, " = ", val)
                else:
                    print("     ", key, " = ", val)
        if (args.caliperConfig):
            print("  *  caliperConfig = ", args.caliperConfig)
        if (args.caliperFilename):
            print("  *  caliperFilename = ", args.caliperFilename)
        if (args.caliperConfigJSON):
            print("  *  caliperConfigJSON = ", args.caliperConfigJSON)
    # Set all the variables.
    gd = globalFrame().f_globals
    for key, val in arg_dict.items():
        if key in options:
            if (type(val) != type(options[key])):
                val = eval(val, gd)
        gd[key] = val
    # Initialize Caliper ConfigManager
    InitTimers(args.caliperConfig,
               args.caliperFilename,
               args.caliperConfigJSON,
               args.caliperOutputDir)
    return

def InitTimers(caliper_config, filename, caliper_json):
    if(caliper_json):
        TimerMgr.load(caliper_json)
        if(not caliper_config):
            raise RuntimeError("SpheralOptionParser: specifying a configuration file without using one of the configurations means no timers are started")
    off_tests = ["none", "off", "disable", "disabled", "0"]
    # Check if Caliper is turned off
    if (caliper_config.lower() in off_tests):
        return
    elif (caliper_config):
        TimerMgr.add(caliper_config)
        TimerMgr.start()
    else:
        import os, sys
        if (filename):
            testname = filename
        else:
            from datetime import datetime
            # Append the current day and time to the filename
            unique_digits = datetime.now().strftime("_%Y_%m_%d_%H%M%S_%f")
            # Name file based on name of python file being run
            testname = os.path.splitext(os.path.basename(sys.argv[0]))[0]
            testname += unique_digits + ".cali"
        TimerMgr.default_start(testname)
    adiak_valueInt("threads_per_rank", omp_get_num_threads())
    adiak_valueInt("num_ranks", mpi.procs)
    return
