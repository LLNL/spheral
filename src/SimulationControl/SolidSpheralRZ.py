#-------------------------------------------------------------------------------
# Import SolidSpheral objects, setting the 2D objects as generic names.
# This version is specialized for the 2D RZ formulation.
#-------------------------------------------------------------------------------
import SolidSpheral
for name in [x for x in SolidSpheral.__dict__ if "2d" in x]:
    exec("%s = SolidSpheral.__dict__['%s']" % (name.replace("2d", ""), name))
for name in [x for x in SolidSpheral.__dict__ if x[-2:] == "RZ"]:
    exec("%s = SolidSpheral.__dict__['%s']" % (name.replace("RZ", ""), name))
del x, name
from SolidSpheral import *
