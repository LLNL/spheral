#-------------------------------------------------------------------------------
# Import Spheral objects, setting the 1-D objects as generic names.
#-------------------------------------------------------------------------------
import Spheral
for name in [x for x in Spheral.__dict__ if "1d" in x]:
    exec("%s = Spheral.__dict__['%s']" % (name.replace("1d", ""), name))
del name
from Spheral import *
FacetedVolume = Spheral.Box1d
