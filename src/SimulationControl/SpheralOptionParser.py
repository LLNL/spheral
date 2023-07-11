#-------------------------------------------------------------------------------
# Create a standard and hopefully convenient command  line parser for Spheral
# scripts.
#-------------------------------------------------------------------------------
import optparse

from SpheralCompiledPackages import *

from SpheralTestUtilities import globalFrame

def commandLine(**options):

    # Build a command line parser with the keyword arguments passed to us.
    parser = optparse.OptionParser()
    for key in options:
        parser.add_option("--" + key,
                          dest = key,
                          default = options[key])

    # Add the universal options suppoted by all Spheral++ scripts.
    parser.add_option("-v", "--verbose",
                      action = "store_true",
                      dest = "__verbose",
                      default = False,
                      help = "Verbose output -- print all options that were set.")

    # Evaluate the command line.
    opts, args = parser.parse_args()

    # Verbose output?
    if opts.__verbose:
        print("All parameters set:")
        for key in options:
            val = eval("opts.%s" % key)
            if val != options[key]:
                print("  *  ", key, " = ", val)
            else:
                print("     ", key, " = ", val)
                
    # Set all the variables.
    gd = globalFrame().f_globals
    for key in options:
        val = eval("opts.%s" % key)

        # The following bit of goofiness is necessary because optparse returns most
        # arguments as a string, but we need the actual types here.
        if type(val) != type(options[key]):
            val = eval(val, gd)

        gd[key] = val

    return
