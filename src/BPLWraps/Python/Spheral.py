# Master import file to import the Spheral packages and a standard set of
# helper extensions.

# Modified version to be compatible with the Boost.Python version of Spheral++.

# ------------------------------------------------------------------------------
# OK, this is nuts.  Check to see if we're built with gcc.  If so, we have to
# do some magic 'cause type info and dynamic_casts get screwed up between
# g++ built libraries.  Sigh.
# ------------------------------------------------------------------------------
try:
    #import sys, dl
    #sys.setdlopenflags(dl.RTLD_NOW|dl.RTLD_GLOBAL)
    import sys
    sys.setdlopenflags(0x100|0x2)
except:
    print "WARNING: unable to set python dl flags on Spheral import."
    pass

# ------------------------------------------------------------------------------
# Import the compiled packages.
# ------------------------------------------------------------------------------
from Geometry import *
from FileIO import *
from DataOutput import *
from NodeList import *
from NodeIterators import *
from Field import *
from FieldOperations import *
from Kernel import *
from SplineKernel import *
from Material import *
from Neighbor import *
from DataBase import *
from Boundary import *
from ArtificialViscosity import *
from Physics import *
from Hydro import *
from ExternalForce import *
from Gravity import *
from Integrator import *
from CXXTypes import *
from Utilities import *

# SPHGravity, baby!
try:
    from SPHGravity import *
except:
    pass

# MHD, too.
try:
    from MHD import *
except:
    pass

# ------------------------------------------------------------------------------
# Import the various FluidNodeLists.
# ------------------------------------------------------------------------------
from FluidNodeLists import *

# ------------------------------------------------------------------------------
# For convenience, shorten the names of the Vector and Tensor classes.
# ------------------------------------------------------------------------------
Vector1d = GeomVector1d
Vector2d = GeomVector2d
Vector3d = GeomVector3d

Tensor1d = GeomTensor1d
Tensor2d = GeomTensor2d
Tensor3d = GeomTensor3d

SymTensor1d = GeomSymTensor1d
SymTensor2d = GeomSymTensor2d
SymTensor3d = GeomSymTensor3d

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
from ExtendFlatFileIO import *

# ------------------------------------------------------------------------------
# If we're running parallel, we need to import the Distributed package as well.
# ------------------------------------------------------------------------------
try:
    import mpi
    mpi.synchronizeQueuedOutput("/dev/null")  # Only have rank 0 print by default.
    try:
        from Distributed import *
    except:
        print 'Warning: you appear to be running pyMPI on ', mpi.procs, \
              'processors, but have not configured Spheral++ for parallel runs.'
        pass
except:
    pass

# ------------------------------------------------------------------------------
# Import the controller and a standard timer class.
# ------------------------------------------------------------------------------
from SpheralTimer import *
from SpheralController import *

# ------------------------------------------------------------------------------
# Initialize the Tau timer package.
# ------------------------------------------------------------------------------
initializeTau()

# ------------------------------------------------------------------------------
# Import the command line interpreter.
# ------------------------------------------------------------------------------
from SpheralOptionParser import commandLine
