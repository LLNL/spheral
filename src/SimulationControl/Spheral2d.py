#-------------------------------------------------------------------------------
# Import Spheral objects, setting the 2-D objects as generic names.
#-------------------------------------------------------------------------------
import Spheral
for name in [x for x in Spheral.__dict__ if "2d" in x]:
    exec("%s = Spheral.__dict__['%s']" % (name.replace("2d", ""), name))
from Spheral import *
FacetedVolume = Spheral.Polygon
