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
    parser.add_argument("--caliper-config", default="", type=str)
    # Evaluate the command line.
    args = parser.parse_args()
    arg_dict = vars(args)

    # Verbose output?
    if args.verbose:
        print("All parameters set:")
        for key, val in arg_dict.items():
            if key in options:
                if val != options[key]:
                    print("  *  ", key, " = ", val)
                else:
                    print("     ", key, " = ", val)
                
    # Set all the variables.
    gd = globalFrame().f_globals
    for key, val in arg_dict.items():
        if key in options:
            if (type(val) != type(options[key])):
                val = eval(val, gd)
            gd[key] = val
    off_tests = ["none", "off", "disable", "disabled"]
    if (not args.caliper_config.lower() in off_tests):
        InitTimers(args)
    return

def InitTimers(args):
    if (args.caliper_config):
        TimerMgr.add(args.caliper_config)
        TimerMgr.start()
    else:
        import random, os, sys
        unique_digits = ''.join(random.sample('0123456789', 4))
        testname = os.path.splitext(os.path.basename(sys.argv[0]))[0] + "_" +  unique_digits
        TimerMgr.default_start(testname)
    return
