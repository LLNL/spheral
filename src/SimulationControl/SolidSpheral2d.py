#-------------------------------------------------------------------------------
# Import SolidSpheral objects, setting the 2-D objects as generic names.
#-------------------------------------------------------------------------------
import SolidSpheral
for name in [x for x in SolidSpheral.__dict__ if "2d" in x]:
    exec("%s = SolidSpheral.__dict__['%s']" % (name.replace("2d", ""), name))
del x, name
from SolidSpheral import *
