# Master import file to import the Spheral packages and a standard set of
# helper extensions.

# Modified version to be compatible with the pybindgen version of Spheral++.

# ------------------------------------------------------------------------------
# OK, this is nuts.  Check to see if we're built with gcc.  If so, we have to
# do some magic 'cause type info and dynamic_casts get screwed up between
# g++ built libraries.  Sigh.
# ------------------------------------------------------------------------------
import sys, ctypes
try:
    sys.setdlopenflags(sys.getdlopenflags() | ctypes.RTLD_GLOBAL)
    #sys.setdlopenflags(ctypes.RTLD_NOW | ctypes.RTLD_GLOBAL)
    #import sys, DLFCN
    #sys.setdlopenflags(sys.getdlopenflags() | DLFCN.RTLD_GLOBAL)
    #sys.setdlopenflags(DLFCN.RTLD_NOW|DLFCN.RTLD_GLOBAL)
    #import sys
    #sys.setdlopenflags(0x100|0x2)
except:
    print "WARNING: unable to set python dl flags on Spheral import."
    pass

# ------------------------------------------------------------------------------
# Import the core Spheral compiled packages.
# ------------------------------------------------------------------------------
from SpheralModules import *
from SpheralModules.Spheral import *
from SpheralModules.Spheral.PythonBoundFunctors import *

# ------------------------------------------------------------------------------
# PolyCliper
# ------------------------------------------------------------------------------
import SpheralModules.PolyClipper as PolyClipper

# ------------------------------------------------------------------------------
# Load up MPI.
# ------------------------------------------------------------------------------
import mpi

# ------------------------------------------------------------------------------
# Import the Material python extensions.
# ------------------------------------------------------------------------------
from MaterialUnits import *
from MaterialEquationsOfState import *

# ------------------------------------------------------------------------------
# Import the various FluidNodeLists.
# ------------------------------------------------------------------------------
from FluidNodeLists import *
from VoidNodeLists import *

# ------------------------------------------------------------------------------
# Import SPH, SVPH, and CRKSPH
# ------------------------------------------------------------------------------
from SPHHydros import *
from PSPHHydros import *
from SVPHHydros import *
from CRKSPHHydros import *
from SPHUtilities import *

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
# If we're running parallel, we need to import the Distributed package as well.
# ------------------------------------------------------------------------------
# try:
#     from SpheralModules.Spheral.PartitionSpace import *
# except:
#     print "Warning: unable to load Distributed components, parallel mode not available."
#     pass

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
    from PolytopeModules import polytope
except:
    print "WARNING: unable to import polytope python bindings."


# ------------------------------------------------------------------------------
# Prepare for timing
# ------------------------------------------------------------------------------
# EasyProfilerStart()

# ------------------------------------------------------------------------------
# Output some useful Spheral configuration info to stdout.
# ------------------------------------------------------------------------------
# boxchars = {'A' : '╔',
#             'B' : '╗',
#             'C' : '╚',
#             'D' : '╝',
#             'E' : '═',
#             'F' : '║'}
# print boxchars['A'] + 78*boxchars['E'] + boxchars['B']
# print "%s  %-76s%s" % (boxchars['F'], "Spheral version: @spheralversion@", boxchars['F'])
# print "%s  %-76s%s" % (boxchars['F'], "  number of MPI tasks       : " + str(mpi.procs), boxchars['F'])
# print "%s  %-76s%s" % (boxchars['F'], "  number of threads per rank: " + str(omp_get_num_threads()), boxchars['F'])
# print boxchars['C'] + 78*boxchars['E'] + boxchars['D']
# print u"╔══════════════════════════════════════════════════════════════════════════════╗"
# print u"║  %-76s║" % "Spheral version: @spheralversion@"
# print u"║  %-76s║" % ("  number of MPI tasks       : " + str(mpi.procs))
# print u"║  %-76s║" % ("  number of threads per rank: " + str(omp_get_num_threads()))
# print u"╚══════════════════════════════════════════════════════════════════════════════╝"
print "/------------------------------------------------------------------------------\\"
print "|  %-76s|" % "Spheral version: @spheralversion@"
print "|  %-76s|" % ("  number of MPI tasks       : " + str(mpi.procs))
print "|  %-76s|" % ("  number of threads per rank: " + str(omp_get_num_threads()))
print "\\------------------------------------------------------------------------------/"
sys.ps1 = "Spheral> "
