"""
Spheral Strength module.

Provides strength modeling helpers.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"Strength/SolidFieldNames.hh"',
                  '"Geometry/Dimension.hh"',
                  '"Geometry/GeomPlane.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# SolidFieldNames
#-------------------------------------------------------------------------------
class SolidFieldNames:

    deviatoricStress = PYB11readonly(static=True)
    deviatoricStressTT = PYB11readonly(static=True)
    plasticStrain = PYB11readonly(static=True)
    plasticStrainRate = PYB11readonly(static=True)
    scalarDamage = PYB11readonly(static=True)
    tensorDamage = PYB11readonly(static=True)
    effectiveTensorDamage = PYB11readonly(static=True)
    damageGradient = PYB11readonly(static=True)
    damageHat = PYB11readonly(static=True)
    strain = PYB11readonly(static=True)
    strainTensor = PYB11readonly(static=True)
    effectiveStrainTensor = PYB11readonly(static=True)
    bulkModulus = PYB11readonly(static=True)
    shearModulus = PYB11readonly(static=True)
    YoungsModulus = PYB11readonly(static=True)
    longitudinalSoundSpeed = PYB11readonly(static=True)
    yieldStrength = PYB11readonly(static=True)
    flaws = PYB11readonly(static=True)
    effectiveFlaws = PYB11readonly(static=True)
    porosityAlpha = PYB11readonly(static=True)
    porosityStrain = PYB11readonly(static=True)
    porosityAlpha0 = PYB11readonly(static=True)
    porosityc0 = PYB11readonly(static=True)
    fragmentIDs = PYB11readonly(static=True)
    particleTypes = PYB11readonly(static=True)
    meltSpecificEnergy = PYB11readonly(static=True)
