#ATS:APGM1 = test(SELF,       " ", np=1, label="Approximate polyhedral gravity model acceleration test")

#-------------------------------------------------------------------------------
# Set up a pair of equal mass N-body points in a simple circular orbit of each
# other.
#-------------------------------------------------------------------------------
from Spheral3d import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from NodeHistory import *
from PolyhedronFileUtilities import *
from SpheralVisitDump import dumpPhysicsState
from math import *
import os, sys

print "Test ApproximatePolyhedralGravityModel"

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(fileName = "data/Bennu_996.obj",  # bennu surface definition
            M = 7.8e10,                        # mass of bennu
            G = 6.67e-11,

            pos = Vector(500.0,0.0,0.0),
            nominalAcceleration = Vector(-0.0000214555067869693,   0.0000000066680564648,  -0.0000000589458328628), 
            errorThreshold = 1e-4,
    )


G = 6.67e-11

# read our mesh
polyBennu = readPolyhedronOBJ(fileName)

# create the gravity model
GravModel = ApproximatePolyhedralGravityModel(polyBennu,M,G)

a = GravModel.acceleration(pos)

error = 100*(a - nominalAcceleration).magnitude()/nominalAcceleration.magnitude()


print "error acceleration field precentage : %0.15f" % error

if error > errorThreshold:
    raise ValueError, "error bounds exceeded"


