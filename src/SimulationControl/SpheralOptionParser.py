#-------------------------------------------------------------------------------
# Create a standard and hopefully convenient command line parser for Spheral
# scripts.
#-------------------------------------------------------------------------------
import argparse, mpi

from SpheralCompiledPackages import *

from SpheralTestUtilities import globalFrame
import SpheralTimingParser

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

    # Parse Caliper and Adiak inputs
    SpheralTimingParser.add_timing_args(parser)

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
    # Initialize timers and add inputs as Adiak metadata
    SpheralTimingParser.init_timer(args)
    return
