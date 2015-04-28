from math import *
from AcousticWaveSolution import AcousticWaveSolution

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
        #self.plus  = AcousticWaveSolution(eos, cs, 0.5*rho0, x0 - 0.5*self.L, x1 + 0.5*self.L,  A,  2.0*pi*k, 0.5*h0)
        #self.minus = AcousticWaveSolution(eos, cs, 0.5*rho0, x0 - 0.5*self.L, x1 + 0.5*self.L, -A, -2.0*pi*k, 0.5*h0)
        return

    # Compute and return the solution on the given positions.
    def solution(self, time, xvals):
        cs2 = self.cs * self.cs
        f1 = -cos(pi*time/(self.L*self.cs))
        f2 = sin(pi*time/(self.L*self.cs))
        v = [self.cs*self.A*sin(pi*self.k*(xi - self.x0)/self.L) * f1 for xi in xvals]
        rho = [self.rho0*(1.0 + self.A*cos(pi*self.k*(xi - self.x0)/self.L) * f2) for xi in xvals]
        u = [cs2*(self.A*cos(pi*self.k*(xi - self.x0)/self.L) * f2) for xi in xvals]
        P = [self.eos.pressure(rhoi, ui) for (rhoi, ui) in zip(rho, u)]
        h = [self.h0*self.rho0/rhoi for rhoi in rho]
        return xvals, v, u, rho, P, h
