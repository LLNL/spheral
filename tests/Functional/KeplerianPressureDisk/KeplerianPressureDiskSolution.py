#-------------------------------------------------------------------------------
# Calcuate the radial profile for the Keplerian disk with arbitrary mix of
# pressure and rotational support.
#-------------------------------------------------------------------------------
from math import *

class KeplerianPressureDiskSolution:

    def __init__(self,
                 fractionPressureSupport,
                 polytropicConstant,
                 polytropicIndex,
                 G,
                 M,
                 coreRadius,
                 rmin = 0.0,
                 rmax = 1.0,
                 nPoints = 101):

        self.fractionPressureSupport = fractionPressureSupport
        self.polytropicConstant = polytropicConstant
        self.polytropicIndex = polytropicIndex
        self.G = G
        self.M = M
        self.coreRadius2 = coreRadius*coreRadius
        self.rmin = rmin
        self.rmax = rmax
        self.nPoints = nPoints

        self.gamma = (polytropicIndex + 1)/polytropicIndex
        return

    def solution(self, time):

        # Return values.
        r = []
        v = []
        u = []
        rho = []
        P = []

        dr = (self.rmax - self.rmin)/(self.nPoints - 1)
        velSupport = 1.0 - self.fractionPressureSupport
        for i in xrange(self.nPoints):
            ri = i*dr
            rhoi = (self.G*self.M/self.polytropicConstant *
                    (self.gamma - 1.0)/self.gamma /
                    sqrt(ri*ri + self.coreRadius2))**(1.0/(self.gamma - 1.0))
            Pi = (self.fractionPressureSupport*self.polytropicConstant*
                  rhoi**self.gamma)
            r.append(ri)
            v.append(ri*sqrt(velSupport*self.G*self.M*(ri*ri + self.coreRadius2)**(-1.5)))
            u.append(Pi/((self.gamma - 1.0)*rhoi))
            rho.append(rhoi)
            P.append(Pi)

        return r, v, u, rho, P
