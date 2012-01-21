#-------------------------------------------------------------------------------
# Import Spheral objects, setting the 3-D objects as generic names.
#-------------------------------------------------------------------------------
import Spheral
for name in [x for x in Spheral.__dict__ if "3d" in x]:
    exec("%s = Spheral.__dict__['%s']" % (name.replace("3d", ""), name))
del x, name
