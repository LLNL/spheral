"""
Spheral Hydro module.

Provides the support classes for hydro algorithms.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"Geometry/GeomPlane.hh"',
                  '"Hydro/HydroFieldNames.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# HydroFieldNames
#-------------------------------------------------------------------------------
class HydroFieldNames:

    mass = PYB11readonly(static=True, returnpolicy="copy")
    position = PYB11readonly(static=True, returnpolicy="copy")
    velocity = PYB11readonly(static=True, returnpolicy="copy")
    H = PYB11readonly(static=True, returnpolicy="copy")
    work = PYB11readonly(static=True, returnpolicy="copy")
    velocityGradient = PYB11readonly(static=True, returnpolicy="copy")
    internalVelocityGradient = PYB11readonly(static=True, returnpolicy="copy")
    hydroAcceleration = PYB11readonly(static=True, returnpolicy="copy")
    massDensity = PYB11readonly(static=True, returnpolicy="copy")
    normalization = PYB11readonly(static=True, returnpolicy="copy")
    specificThermalEnergy = PYB11readonly(static=True, returnpolicy="copy")
    maxViscousPressure = PYB11readonly(static=True, returnpolicy="copy")
    effectiveViscousPressure = PYB11readonly(static=True, returnpolicy="copy")
    massDensityCorrection = PYB11readonly(static=True, returnpolicy="copy")
    viscousWork = PYB11readonly(static=True, returnpolicy="copy")
    XSPHDeltaV = PYB11readonly(static=True, returnpolicy="copy")
    XSPHWeightSum = PYB11readonly(static=True, returnpolicy="copy")
    Hsmooth = PYB11readonly(static=True, returnpolicy="copy")
    massFirstMoment = PYB11readonly(static=True, returnpolicy="copy")
    massSecondMoment = PYB11readonly(static=True, returnpolicy="copy")
    weightedNeighborSum = PYB11readonly(static=True, returnpolicy="copy")
    pressure = PYB11readonly(static=True, returnpolicy="copy")
    temperature = PYB11readonly(static=True, returnpolicy="copy")
    soundSpeed = PYB11readonly(static=True, returnpolicy="copy")
    pairAccelerations = PYB11readonly(static=True, returnpolicy="copy")
    pairWork = PYB11readonly(static=True, returnpolicy="copy")
    gamma = PYB11readonly(static=True, returnpolicy="copy")
    entropy = PYB11readonly(static=True, returnpolicy="copy")
    PSPHcorrection = PYB11readonly(static=True, returnpolicy="copy")
    omegaGradh = PYB11readonly(static=True, returnpolicy="copy")
    numberDensitySum = PYB11readonly(static=True, returnpolicy="copy")
    timeStepMask = PYB11readonly(static=True, returnpolicy="copy")
    surfacePoint = PYB11readonly(static=True, returnpolicy="copy")
    voidPoint = PYB11readonly(static=True, returnpolicy="copy")
    etaVoidPoints = PYB11readonly(static=True, returnpolicy="copy")
    cells = PYB11readonly(static=True, returnpolicy="copy")
    cellFaceFlags = PYB11readonly(static=True, returnpolicy="copy")
    M_SPHCorrection = PYB11readonly(static=True, returnpolicy="copy")
    volume = PYB11readonly(static=True, returnpolicy="copy")
    linearMomentum = PYB11readonly(static=True, returnpolicy="copy")
    totalEnergy = PYB11readonly(static=True, returnpolicy="copy")
    mesh = PYB11readonly(static=True, returnpolicy="copy")
    hourglassMask = PYB11readonly(static=True, returnpolicy="copy")
    faceVelocity = PYB11readonly(static=True, returnpolicy="copy")
    faceForce = PYB11readonly(static=True, returnpolicy="copy")
    faceMass = PYB11readonly(static=True, returnpolicy="copy")
    polyvols = PYB11readonly(static=True, returnpolicy="copy")
    massDensityGradient = PYB11readonly(static=True, returnpolicy="copy")
    specificHeat = PYB11readonly(static=True, returnpolicy="copy")
    normal = PYB11readonly(static=True, returnpolicy="copy")
    surfaceArea = PYB11readonly(static=True, returnpolicy="copy")
