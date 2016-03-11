#-------------------------------------------------------------------------------
# Import SolidSpheral objects, setting the 1-D objects as generic names.
#-------------------------------------------------------------------------------
import SolidSpheral
for name in [x for x in SolidSpheral.__dict__ if "1d" in x]:
    exec("%s = SolidSpheral.__dict__['%s']" % (name.replace("1d", ""), name))
del x, name
from SolidSpheral import *
