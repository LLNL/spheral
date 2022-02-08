from SpheralCompiledPackages import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic SPHHydro pattern.
#-------------------------------------------------------------------------------
DEMFactoryString = """
class %(classname)s%(dim)s(DEMBase%(dim)s):

    def __init__(self,
                 dataBase,
                 W,
                 stepsPerCollision = 25,
                 xmin = Vector%(dim)s(-1e100, -1e100, -1e100),
                 xmax = Vector%(dim)s( 1e100,  1e100,  1e100)):

        DEMBase%(dim)s.__init__(self,
                                dataBase,
                                W,
                                stepsPerCollision,
                                xmin,
                                xmax)
        return
"""

#-------------------------------------------------------------------------------
# Make 'em.
#-------------------------------------------------------------------------------
for dim in dims:
    exec(DEMFactoryString % {"dim"                  : "%id" % dim,
                             "classname"            : "DEM"})

#-------------------------------------------------------------------------------
# Provide a factory function to return the appropriate SPH hydro.
#-------------------------------------------------------------------------------
def DEM(dataBase,
        W = None,
        stepsPerCollision = 25,
        xmin = (-1e100, -1e100, -1e100),
        xmax = ( 1e100,  1e100,  1e100)):

    assert dataBase.numDEMNodeLists == dataBase.numNodeLists, "all nodelists must be dem nodelists"
    assert stepsPerCollision > 10, "stepsPerCollision too low, reccomended is 25-50"

    ndim = dataBase.nDim

    Constructor = eval("DEM%id" % ndim)

    # use a course c2 kernel by default (this is just for the neighbor search)
    if W is None:
        W = eval("TableKernel%id(WendlandC2Kernel%id(), 10)" % (ndim,ndim))

    # Build the constructor arguments
    xmin = (ndim,) + xmin
    xmax = (ndim,) + xmax
    kwargs = {"dataBase" : dataBase,
              "W" : W,
              "stepsPerCollision" : stepsPerCollision,
              "xmin" : eval("Vector%id(%g, %g, %g)" % xmin),
              "xmax" : eval("Vector%id(%g, %g, %g)" % xmax)}

    # Build and return the thing.
    result = Constructor(**kwargs)

    return result




