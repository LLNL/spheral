from SpheralCompiledPackages import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic SPHHydro pattern.
#-------------------------------------------------------------------------------
DEMFactoryString = """
class %(classname)s%(dim)s(LinearSpringDEM%(dim)s):

    def __init__(self,
                 dataBase,
                 normalSpringConstant,
                 restitutionCoefficient,
                 stepsPerCollision,
                 xmin = Vector%(dim)s(-1e100, -1e100, -1e100),
                 xmax = Vector%(dim)s( 1e100,  1e100,  1e100)):

        LinearSpringDEM%(dim)s.__init__(self,
                                        dataBase,
                                        normalSpringConstant,
                                        restitutionCoefficient,
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
def LinearSpringDEM(dataBase,
                    normalSpringConstant,
                    restitutionCoefficient,
                    stepsPerCollision = 25,
                    xmin = (-1e100, -1e100, -1e100),
                    xmax = ( 1e100,  1e100,  1e100)):
    assert dataBase.numDEMNodeLists == dataBase.numNodeLists, "all nodelists must be dem nodelists"
    assert stepsPerCollision > 10, "stepsPerCollision too low, reccomended is 25-50"

    ndim = dataBase.nDim

    Constructor = eval("LinearSpringDEM%id" % ndim)

    # Build the constructor arguments
    xmin = (ndim,) + xmin
    xmax = (ndim,) + xmax
    kwargs = {"dataBase" : dataBase,
              "normalSpringConstant" : normalSpringConstant,
              "restitutionCoefficient" : restitutionCoefficient,
              "stepsPerCollision" : stepsPerCollision,
              "xmin" : eval("Vector%id(%g, %g, %g)" % xmin),
              "xmax" : eval("Vector%id(%g, %g, %g)" % xmax)}

    # Build and return the thing.
    result = Constructor(**kwargs)

    return result

#-------------------------------------------------------------------------------
# convienence function that defaults to Linear Spring DEM
#-------------------------------------------------------------------------------
def DEM(dataBase,
        normalSpringConstant,
        restitutionCoefficient,
        stepsPerCollision = 25,
        xmin = (-1e100, -1e100, -1e100),
        xmax = ( 1e100,  1e100,  1e100)):
    return LinearSpringDEM(dataBase,
                           normalSpringConstant,
                           restitutionCoefficient,
                           stepsPerCollision,
                           xmin,
                           xmax)




