from math import *
from SpheralTestUtilities import *

#-------------------------------------------------------------------------------
# A method to adjust the node positions and masses of the cylinder to match the
# beveling used in the experiments.
#
# A subtly here is that in 2-D we assume that the tube is aligned along the x
# axis, but in 3-D along the z.  :(
#-------------------------------------------------------------------------------
def bevelTubeEntrance(tubeNodeList,
                      nDim,            # number of dimensions (2,3)
                      openingAngle,    # radians
                      tubeInnerRadius, # length units
                      tubeThickness,   # length units
                      xBevelBegin):    # length units

    # Pre-conditions.
    assert nDim == 2 or nDim == 3
    assert openingAngle >= 0.0
    assert tubeInnerRadius >= 0.0
    assert tubeThickness > 0.0

    # Pre-compute some useful geometry.
    tantheta = tan(openingAngle)
    tubeOuterRadius = tubeInnerRadius + tubeThickness
    vol0 = -1.0
    if nDim == 2:
        vol0 = tubeThickness
    else:
        vol0 = pi*(tubeOuterRadius**2 - tubeInnerRadius**2)
    assert vol0 > 0.0

    # Get the position and mass fields.
    position = tubeNodeList.positions()
    mass = tubeNodeList.mass()

    # Iterate over the nodes in the node list.
    numNodesBeveled = 0
    for i in xrange(tubeNodeList.numInternalNodes):

        # Is this node in the x-range to be beveled?
        xi = position[i].x
        yi = position[i].y
        zi = 0.0
        if nDim == 3:
            xi = position[i].z
            zi = position[i].x
        if xi > xBevelBegin:
            numNodesBeveled += 1

            # Adjust the position.
            dThickness = tantheta * (xi - xBevelBegin)
            assert dThickness >= 0.0 and dThickness < tubeThickness
            rmult = 1.0 - dThickness/tubeThickness
            assert rmult > 0.0 and rmult <= 1.0
            if nDim == 2:
                assert distinctlyGreaterThan(yi, 0.0)
                drFromOuterRadius = rmult*(tubeOuterRadius - yi)
                assert drFromOuterRadius >= 0.0 and drFromOuterRadius <= tubeOuterRadius
                dr = tubeOuterRadius - drFromOuterRadius
                A = dr/yi
                assert A >= 1.0
                position[i].y *= A
                assert position[i].y >= yi and position[i].y <= tubeOuterRadius
            else:
                drold = sqrt(yi*yi + zi*zi)
                assert distinctlyGreaterThan(drold, 0.0)
                drFromOuterRadius = rmult*(tubeOuterRadius - drold)
                assert drFromOuterRadius >= 0.0 and drFromOuterRadius <= tubeOuterRadius
                dr = tubeOuterRadius - drFromOuterRadius
                A = dr/drold
                assert A >= 1.0
                position[i].x *= A
                position[i].y *= A
                drnew = sqrt(position[i].x**2 + position[i].y**2)
                assert drnew >= drold and drnew <= tubeOuterRadius

            # Adjust the node mass.
            dvol = -1.0
            if nDim == 2:
                dvol = dThickness
            else:
                dvol = pi*((tubeInnerRadius + dThickness)**2 - tubeInnerRadius**2)
            assert dvol >= 0.0
            dm = dvol/vol0
            mmult = 1.0 - dm
            mass[i] *= mmult

    return numNodesBeveled

