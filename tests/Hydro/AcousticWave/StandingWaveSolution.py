from math import *
from Spheral1d import *

#-------------------------------------------------------------------------------
# The analytic answer for the standing wave problem.
#-------------------------------------------------------------------------------
class StandingWaveSolution:

    # Constructor.
    def __init__(self,
                 eos,
                 cs,
                 rho0,
                 x0,
                 x1,
                 A,
                 k,
                 h0):
        self.eos = eos
        self.cs = cs
        self.rho0 = rho0
        self.x0 = x0
        self.x1 = x1
        self.A = A
        self.k = k
        self.h0 = h0
        self.L = x1 - x0
        return

    # Compute and return the solution on the given positions.
    def solution(self, time, xvals):
        cs2 = self.cs * self.cs
        f1 = -cos(pi*time/(self.L*self.cs))
        f2 = sin(pi*time/(self.L*self.cs))
        v = [self.cs*self.A*sin(pi*self.k*(xi - self.x0)/self.L) * f1 for xi in xvals]
        rho = [self.rho0*(1.0 + self.A*cos(pi*self.k*(xi - self.x0)/self.L) * f2) for xi in xvals]
        u = [cs2*(self.A*cos(pi*self.k*(xi - self.x0)/self.L) * f2) for xi in xvals]
        h = [self.h0*self.rho0/rhoi for rhoi in rho]

        n = len(xvals)
        nodes = makeFluidNodeList("tmp nodes", self.eos, 
                                  numInternal = n)
        rhof = ScalarField("density", nodes)
        uf = ScalarField("sp energy", nodes)
        for i in xrange(n):
            rhof[i] = rho[i]
            uf[i] = u[i]
        Pf = ScalarField("pressure", nodes)
        self.eos.setPressure(Pf, rhof, uf)
        P = list(Pf.internalValues())

        return xvals, v, u, rho, P, h
