# Master import file to import the Spheral packages and a standard set of
# helper extensions.

# Modified version to be compatible with the pybindgen version of Spheral++.

from SpheralUtilities import BuildData

if not BuildData.cxx_compiler_id == "GNU":
    try:
        import sys, ctypes
        sys.setdlopenflags(sys.getdlopenflags() | ctypes.RTLD_GLOBAL)
    except:
        print "WARNING: unable to set python dl flags on Spheral import."
        pass

# ------------------------------------------------------------------------------
# Load up MPI.
# ------------------------------------------------------------------------------
import mpi

# ------------------------------------------------------------------------------
# Import the compiled packages.
# ------------------------------------------------------------------------------
from SpheralCompiledPackages import *

# ------------------------------------------------------------------------------
# Import the Material python extensions.
# ------------------------------------------------------------------------------
from MaterialUnits import *
from MaterialEquationsOfState import *

# ------------------------------------------------------------------------------
# Import the various NodeLists.
# ------------------------------------------------------------------------------
from FluidNodeLists import *
from SolidNodeLists import *
from VoidNodeLists import *
from DEMNodeLists import *

# ------------------------------------------------------------------------------
# Import SPH, SVPH, and CRKSPH
# ------------------------------------------------------------------------------
from SPHHydros import *
from PSPHHydros import *
from GSPHHydros import *
from FSISPHHydros import *
from SlideSurfaces import *
#from SVPHHydros import *
from CRKSPHHydros import *
#from TaylorSPHHydros import *
from DEM import *
from SPHUtilities import *

# ------------------------------------------------------------------------------
# Import the SolidMaterial python extensions.
# ------------------------------------------------------------------------------
from SolidMaterialUnits import *
from SolidMaterialEquationsOfState import *

from GradyKippTensorDamage import *
from JohnsonCookDamageFactories import (JohnsonCookDamageConstant,
                                        JohnsonCookDamageGaussian,
                                        JohnsonCookDamageWeibull)

# ------------------------------------------------------------------------------
# Helpful things with strings.
# ------------------------------------------------------------------------------
## import CXXTypes
## def add_string_str(self, pystr):
##     return str(self).join(pystr)

## def add_str_string(self, cppstr):
##     return self.join(str(cppstr))

## CXXTypes.string.__add__ = add_string_str

# ------------------------------------------------------------------------------
# Import the picking helpers that register how to pickle Spheral++ objects.
# ------------------------------------------------------------------------------
from SpheralPickle import *

# ------------------------------------------------------------------------------
# Import the python enhanced FlatFileIO object.
# ------------------------------------------------------------------------------
#from ExtendFlatFileIO import *

# ------------------------------------------------------------------------------
# Import the controller and a standard timer class.
# ------------------------------------------------------------------------------
from SpheralTimer import *
from SpheralController import *

# ------------------------------------------------------------------------------
# Import the command line interpreter.
# ------------------------------------------------------------------------------
from SpheralOptionParser import commandLine

# ------------------------------------------------------------------------------
# See if we can import the polytope bindings.
# ------------------------------------------------------------------------------
try:
    import polytope
except:
    print "WARNING: unable to import polytope python bindings."

# ------------------------------------------------------------------------------
# Import our shadow layers for augmenting C++ types.
# ------------------------------------------------------------------------------
for shadowedthing in ("TillotsonEquationOfState",
                      "GruneisenEquationOfState",
                      "ConstantStrength",
                      "ProbabilisticDamageModel",
                      "IvanoviSALEDamageModel",
                      "ANEOS"):
    for dim in dims:
        exec("from Shadow%(thing)s import %(thing)s%(dim)sd" % {"thing" : shadowedthing,
                                                                "dim"   : dim})

# ------------------------------------------------------------------------------
# Prepare for timing
# ------------------------------------------------------------------------------
# EasyProfilerStart()

# ------------------------------------------------------------------------------
# Output some useful Spheral configuration info to stdout
# ------------------------------------------------------------------------------
print "/------------------------------------------------------------------------------\\"
import Spheral_banner
print "|  %-76s|" % ("  number of MPI tasks       : " + str(mpi.procs))
print "|  %-76s|" % ("  number of threads per rank: " + str(omp_get_num_threads()))
print "\\------------------------------------------------------------------------------/"

# ------------------------------------------------------------------------------
# Set the prompt just to clear to folks they now have Spheral
# ------------------------------------------------------------------------------
sys.ps1 = "Spheral> "
