#-------------------------------------------------------------------------------
# Create a standard and hopefully convenient command line parser for Spheral
# scripts.
#-------------------------------------------------------------------------------
import argparse

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
    parser.add_argument("--caliperConfig", default="", type=str)
    parser.add_argument("--caliperFilename", default="", type=str)
    # Evaluate the command line.
    args = parser.parse_args()
    arg_dict = vars(args)

    if (not TimerMgr.timers_usable()):
        if (args.caliperConfig or args.caliperFilename):
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
    # Set all the variables.
    gd = globalFrame().f_globals
    for key, val in arg_dict.items():
        if key in options:
            if (type(val) != type(options[key])):
                val = eval(val, gd)
        gd[key] = val
    InitTimers(args.caliperConfig, args.caliperFilename)
    return

def InitTimers(caliper_config, filename):
    off_tests = ["none", "off", "disable", "disabled", "0"]
    if (caliper_config.lower() in off_tests):
        return
    elif (caliper_config):
        TimerMgr.add(caliper_config)
        TimerMgr.start()
    else:
        import random, os, sys
        if (filename):
            testname = filename
            # Remove the cali file extension
            if (".cali" in testname):
                testname = testname.replace(".cali", "")
        else:
            unique_digits = ''.join(random.sample('0123456789', 4))
            testname = os.path.splitext(os.path.basename(sys.argv[0]))[0] + "_" +  unique_digits
        TimerMgr.default_start(testname)
    return
