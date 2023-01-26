from math import *

from NodeGeneratorBase import *
from Spheral import Vector3d, Tensor3d, SymTensor3d, Plane3d, rotationMatrix
from SpheralTestUtilities import *

import mpi
procID = mpi.rank
nProcs = mpi.procs

#-------------------------------------------------------------------------------
# The sign function
#-------------------------------------------------------------------------------
def sgn(x):
    if x >= 0.0:
        return 1.0
    else:
        return -1.0

#-------------------------------------------------------------------------------
# Produce spherical shell sections, bounded by four planes.
#-------------------------------------------------------------------------------
class GenerateSphericalShellSection(NodeGeneratorBase):

    #-------------------------------------------------------------------------------
    # Constructor
    #-------------------------------------------------------------------------------
    def __init__(self,
                 nr,
                 nl,
                 rho0,
                 r0,
                 rthick,
                 openingAngle,
                 nNodePerh = 2.01,
                 SPH = False):
        self.nr = nr
        self.nl = nl
        self.rho0 = rho0
        self.r0 = r0
        self.rthick = rthick
        self.openingAngle = openingAngle
        self.nNodePerh = nNodePerh

        # Prepare our internal arrays.
        self.x = []
        self.y = []
        self.z = []
        self.m = []
        self.H = []

        # Useful stuff to precompute.
        phi0 = -0.5*openingAngle
        dphi = openingAngle/nl
        dr = rthick/nr
        halfpi = 0.5*pi

        # We pretend that the subvolumes per point are pretty much the same,
        # just rotated into a different orientation.  Not a good approximation
        # as we increase the opening angle, but OK for small angles.
        dx = r0*dphi
        V0 = dx*dx*dr
        m0 = V0*rho0
        H0 = SymTensor3d(1.0/dx, 0.0, 0.0,
                         0.0, 1.0/dx, 0.0,
                         0.0, 0.0, 1.0/dr)/nNodePerh

        # Iterate over the angles.
        for ip1 in range(nl):
            phi1 = phi0 + (ip1 + 0.5)*dphi
            plane1 = self.computePlane(0.0, phi1)
            for ip2 in range(nl):
                phi2 = phi0 + (ip2 + 0.5)*dphi
                plane2 = self.computePlane(halfpi, phi2)

                # Get the line that represents the intersection of the two planes.
                P, N = self.computePlaneIntersection(plane1, plane2)

                # Intersect this line with each radii we want, which gives us
                # the desired positions.
                for ir in range(nr):
                    r = r0 + (ir + 0.5)*dr
                    pos = self.intersectLineAndSphere(P, N, r)
                    assert fuzzyEqual(pos.magnitude(), r, 1.0e-10)
                    self.x.append(pos.x)
                    self.y.append(pos.y)
                    self.z.append(pos.z)
                    self.m.append(m0)

                    # Rotate H into the appropriate frame.
                    R = rotationMatrix(N)
                    Hi = SymTensor3d(H0)
                    Hi.rotationalTransform(R)
                    self.H.append(Hi)

        # Check that the points are within bounds.
        planes = self.boundaryPlanes()
        for x, y, z in zip(self.x, self.y, self.z):
            for plane in planes:
                assert Vector3d(x, y, z) > plane

        # Have the base class break up the serial node distribution
        # for parallel cases.
        NodeGeneratorBase.__init__(self, True,
                                   self.x, self.y, self.z, self.m, self.H)

        # If SPH has been specified, make sure the H tensors are round.
        if SPH:
            self.makeHround()

        return

    #-------------------------------------------------------------------------------
    # Return the plane corresponding to the given displacment in theta and phi.
    #-------------------------------------------------------------------------------
    def computePlane(self, theta, phi):
        ct = cos(theta)
        st = sin(theta)
        cp = cos(phi)
        sp = sin(phi)
        P = Vector3d(ct*sp, st*sp, cp)
        N = Vector3d(-ct*cp, -st*cp, sp)
        assert fuzzyEqual(N.magnitude2(), 1.0)
        return Plane3d(P, N)

    #-------------------------------------------------------------------------------
    # Return the line resulting from the intersection of two planes.
    #-------------------------------------------------------------------------------
    def computePlaneIntersection(self, plane1, plane2):
        p1 = plane1.point()
        n1 = plane1.normal()
        p2 = plane2.point()
        n2 = plane2.normal()
        assert distinctlyLessThan(abs(n1.dot(n2)), 1.0)
        d1 = -n1.dot(p1)
        d2 = -n2.dot(p2)
        n3 = (n1.cross(n2)).unitVector()
        p3 = (d2*n1 - d1*n2).cross(n3) / (n1.cross(n2).dot(n3))
        return p3, n3

    #-------------------------------------------------------------------------------
    # Return the positive intersection of a line with a sphere centered at the
    # origin.
    #-------------------------------------------------------------------------------
    def intersectLineAndSphere(self, point, normal, r):
        b = 2.0*normal.dot(point)
        c = point.magnitude2() - r*r
        q = -0.5*(b + sgn(b)*sqrt(b*b - 4.0*c))
        assert distinctlyGreaterThan(abs(q), 0.0)
        a1 = q
        a2 = c/q
        p1 = (point + a1*normal).unitVector() * r
        p2 = (point + a2*normal).unitVector() * r
        if p1 > p2:
            return p1
        else:
            return p2

    #-------------------------------------------------------------------------------
    # Return the boundary planes.
    #-------------------------------------------------------------------------------
    def boundaryPlanes(self):
        p = 0.5*self.openingAngle
        halfpi = 0.5*pi
        return (self.computePlane(0.0, p),
                self.computePlane(halfpi, p),
                self.computePlane(pi, p),
                self.computePlane(-halfpi, p))

    #-------------------------------------------------------------------------------
    # Get the position for the given node index.
    #-------------------------------------------------------------------------------
    def localPosition(self, i):
        assert i >= 0 and i < len(self.x)
        assert len(self.x) == len(self.y)
        return Vector3d(self.x[i], self.y[i], self.z[i])

    #-------------------------------------------------------------------------------
    # Get the mass for the given node index.
    #-------------------------------------------------------------------------------
    def localMass(self, i):
        assert i >= 0 and i < len(self.m)
        return self.m[i]

    #-------------------------------------------------------------------------------
    # Get the mass density for the given node index.
    #-------------------------------------------------------------------------------
    def localMassDensity(self, i):
        assert i >= 0 and i < len(self.x)
        return self.rho0

    #-------------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #-------------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]

