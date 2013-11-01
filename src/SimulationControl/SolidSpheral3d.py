#-------------------------------------------------------------------------------
# Import SolidSpheral objects, setting the 3-D objects as generic names.
#-------------------------------------------------------------------------------
import SolidSpheral
for name in [x for x in SolidSpheral.__dict__ if x[-2:] == "3d"]:
    exec("%s = SolidSpheral.__dict__['%s']" % (name[:-2], name))
del x, name
from SolidSpheral import *
