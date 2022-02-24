#-------------------------------------------------------------------------------
# Import Spheral objects, setting the 1D objects as generic names.
# This version is specialized for the 1D spherical formulation.
#-------------------------------------------------------------------------------
import Spheral
for name in [x for x in Spheral.__dict__ if "1d" in x]:
    exec("%s = Spheral.__dict__['%s']" % (name.replace("1d", ""), name))
for name in [x for x in Spheral.__dict__ if (x != "Spherical" and x[:9] == "Spherical")]:
    exec("%s = Spheral.__dict__['%s']" % (name.replace("Spherical", ""), name))
del x, name
from Spheral import *
GeometryRegistrar.coords(CoordinateType.Spherical)
