# Master import file to import the Spheral packages and a standard set of
# helper extensions.

# Modified version to be compatible with the pybindgen version of Spheral++.

# ------------------------------------------------------------------------------
# OK, this is nuts.  Check to see if we're built with gcc.  If so, we have to
# do some magic 'cause type info and dynamic_casts get screwed up between
# g++ built libraries.  Sigh.
# ------------------------------------------------------------------------------
try:
    import sys, ctypes
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
# Load up MPI.
# ------------------------------------------------------------------------------
import mpi

# ------------------------------------------------------------------------------
# Import the compiled packages.
# ------------------------------------------------------------------------------
from SpheralModules import *
from SpheralModules.Spheral import *
from SpheralModules.Spheral.NodeSpace import *
from SpheralModules.Spheral.FieldSpace import *
from SpheralModules.Spheral.DataBaseSpace import *
from SpheralModules.Spheral.FileIOSpace import *
from SpheralModules.Spheral.ArtificialViscositySpace import *
from SpheralModules.Spheral.DataOutput import *
from SpheralModules.Spheral.KernelSpace import *
from SpheralModules.Spheral.NeighborSpace import *
from SpheralModules.Spheral.Material import *
from SpheralModules.Spheral.BoundarySpace import *
from SpheralModules.Spheral.PhysicsSpace import *
from SpheralModules.Spheral.GravitySpace import *
from SpheralModules.Spheral.IntegratorSpace import *
from SpheralModules.Spheral.SPHSpace import *
from SpheralModules.Spheral.CRKSPHSpace import *
from SpheralModules.Spheral.SVPHSpace import *
from SpheralModules.Spheral.TaylorSPHSpace import *
from SpheralModules.Spheral.MeshSpace import *

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
from SVPHHydros import *
from CRKSPHHydros import *
from TaylorSPHHydros import *
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
try:
    from SpheralModules.Spheral.PartitionSpace import *
except:
    print "Warning: unable to load Distributed components, parallel mode not available."
    pass

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

