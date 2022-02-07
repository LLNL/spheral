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
                 cfl = 0.25,
                 xmin = Vector%(dim)s(-1e100, -1e100, -1e100),
                 xmax = Vector%(dim)s( 1e100,  1e100,  1e100)):

        DEMBase%(dim)s.__init__(self,
                                dataBase,
                                W,
                                cfl,
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
        cfl = 0.05,
        xmin = (-1e100, -1e100, -1e100),
        xmax = ( 1e100,  1e100,  1e100)):

    assert dataBase.numDEMNodeLists == dataBase.numNodeLists

    ndim = dataBase.nDim

    Constructor = eval("DEM%id" % ndim)

    # use a course c2 kernel by default (this is just for the neighbor search)
    if W is None:
        W = eval("TableKernel%id(WendlandC2Kernel2%id(), 10)" % (ndim,ndim))

    # Build the constructor arguments
    xmin = (ndim,) + xmin
    xmax = (ndim,) + xmax
    kwargs = {"dataBase" : dataBase,
              "W" : W,
              "cfl" : cfl,
              "xmin" : eval("Vector%id(%g, %g, %g)" % xmin),
              "xmax" : eval("Vector%id(%g, %g, %g)" % xmax)}

    # Build and return the thing.
    result = Constructor(**kwargs)

    return result




