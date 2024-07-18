#-------------------------------------------------------------------------------
# Adiak utilities
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11cppname("spheral_adiak_init")
def adiak_init():
    "Initialize Adiak and run collect_all"
    return "void"

@PYB11cppname("adiak::fini")
def adiak_fini():
    "Finalize Adiak"
    return "void"
