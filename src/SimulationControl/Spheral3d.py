#-------------------------------------------------------------------------------
# Import Spheral objects, setting the 3-D objects as generic names.
#-------------------------------------------------------------------------------
import Spheral
for name in [x for x in Spheral.__dict__ if x[-2:] == "3d"]:
    exec("%s = Spheral.__dict__['%s']" % (name[:-2], name))
del x, name
from Spheral import *
FacetedVolume = Spheral.Polyhedron
