from SpheralCompiledPackages import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# linear spring helper class
#-------------------------------------------------------------------------------
def LinearSpringDEM(dataBase,
                    normalSpringConstant,
                    normalRestitutionCoefficient,
                    tangentialSpringConstant,
                    tangentialRestitutionCoefficient,
                    dynamicFrictionCoefficient,
                    staticFrictionCoefficient,
                    rollingFrictionCoefficient,
                    torsionalFrictionCoefficient,
                    cohesiveTensileStrength=0.0,
                    shapeFactor=0.0,
                    stepsPerCollision = 25,
                    xmin = (-1e100, -1e100, -1e100),
                    xmax = ( 1e100,  1e100,  1e100)):
    assert dataBase.numDEMNodeLists == dataBase.numNodeLists, "all nodelists must be dem nodelists"
    assert stepsPerCollision > 1, "stepsPerCollision too low, reccomended is 25-50"

    ndim = dataBase.nDim

    Constructor = eval("LinearSpringDEM%id" % ndim)

    # Build the constructor arguments
    xmin = (ndim,) + xmin
    xmax = (ndim,) + xmax
    kwargs = {"dataBase" : dataBase,
              "normalSpringConstant" : normalSpringConstant,
              "normalRestitutionCoefficient" : normalRestitutionCoefficient,
              "tangentialSpringConstant" : tangentialSpringConstant,
              "tangentialRestitutionCoefficient" : tangentialRestitutionCoefficient,
              "dynamicFrictionCoefficient" : dynamicFrictionCoefficient,
              "staticFrictionCoefficient" : staticFrictionCoefficient,
              "rollingFrictionCoefficient" : rollingFrictionCoefficient,
              "torsionalFrictionCoefficient" : torsionalFrictionCoefficient,
              "cohesiveTensileStrength" : cohesiveTensileStrength,
              "shapeFactor" : shapeFactor,
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
        normalRestitutionCoefficient,
        tangentialSpringConstant,
        tangentialRestitutionCoefficient,
        dynamicFrictionCoefficient,
        staticFrictionCoefficient,
        rollingFrictionCoefficient,
        torsionalFrictionCoefficient,
        cohesiveTensileStrength,
        shapeFactor,
        stepsPerCollision = 25,
        xmin = (-1e100, -1e100, -1e100),
        xmax = ( 1e100,  1e100,  1e100)):
    return LinearSpringDEM(dataBase,
                           normalSpringConstant,
                           normalRestitutionCoefficient,
                           tangentialSpringConstant,
                           tangentialRestitutionCoefficient,
                           dynamicFrictionCoefficient,
                           staticFrictionCoefficient,
                           rollingFrictionCoefficient,
                           torsionalFrictionCoefficient,
                           cohesiveTensileStrength,
                           shapeFactor,
                           stepsPerCollision,
                           xmin,
                           xmax)




