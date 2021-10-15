from math import *

#-------------------------------------------------------------------------------
# The analytic answer for the acoustic wave problem.
#-------------------------------------------------------------------------------
class AcousticWaveSolutionMod:

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
        return

    # Compute and return the solution on the given positions.
    def solution(self, time, xvals):
        omegat = self.k*self.cs*time
        length = self.x1 - self.x0
        v = [self.cs*self.A*sin(self.k*(x - self.x0)/length - omegat) for x in xvals]
        u = [1.0]*len(xvals)
        rho = [self.rho0*(1.0 + self.A*sin(self.k*(x - self.x0)/length - omegat)) for x in xvals]
        P = [self.eos.pressure(rhoi, ui) for (rhoi, ui) in zip(rho, u)]
        h = [self.h0*rhoi/self.rho0 for rhoi in rho]
        return xvals, v, u, rho, P, h
