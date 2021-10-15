#-------------------------------------------------------------------------------
# This class provides the self-similar analytic solution to the
# Kidder isentropic ICF capsule impolosion.  See references:
#
# R.E. Kidder, Nucl. Fusion 16 (1976), 3-14.
# Maire, JCP 228 (2009), 6882-6915.
#-------------------------------------------------------------------------------
from math import *
from SpheralTestUtilities import sgn
from numericalIntegration import trapezoidalIntegration

class KidderIsentropicCapsuleAnalyticSolution:

    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self,
                 nu,      # Dimension (1, 2, or 3).
                 r0,      # Initial inner boundary radius
                 r1,      # Initial outer boundary radius
                 P0,      # Initial inner boundary pressure
                 P1,      # Initial outer boundary pressure
                 rho1,    # Initial outer boundary density
                 ):

        # Pre-conditions.
        assert nu in (1, 2, 3)
        assert r0 < r1
        assert rho1 > 0.0

        # Derive a bunch of useful parameters.
        self.nu = nu
        self.r0 = r0
        self.r1 = r1
        self.P0 = P0
        self.P1 = P1
        self.rho1 = rho1
        self.gamma = 1.0 + 2.0/nu
        self.gamma1 = self.gamma - 1.0
        self.S = P1/(rho1**self.gamma)
        self.rho0 = rho1*(P0/P1)**(1.0/self.gamma)
        self.tau = sqrt(0.5*self.gamma1 * (r1*r1 - r0*r0)/(self.S*self.gamma*(self.rho1**self.gamma1 - self.rho0**self.gamma1)))
        assert self.tau > 0.0

        # Time 0 energy constants.
        self.K0 = trapezoidalIntegration(self._K0profile, r0, r1, 2001)
        self.T0 = trapezoidalIntegration(self._T0profile, r0, r1, 2001)

        return

    #---------------------------------------------------------------------------
    # tfrac : The dimensionless time.
    #---------------------------------------------------------------------------
    def tfrac(self, t):
        t1 = t/self.tau
        assert t1 >= 0.0 and t1 <= 1.0
        return t1
        
    #---------------------------------------------------------------------------
    # hfrac : the dimensionless radius scaling.
    #---------------------------------------------------------------------------
    def hfrac(self, t):
        return sqrt(max(0.0, 1.0 - self.tfrac(t)**2))

    #---------------------------------------------------------------------------
    # hfracDot : the time derivative of the dimensionless radius scaling.
    #---------------------------------------------------------------------------
    def hfracDot(self, t):
        return -t/(self.tau*self.tau * self.hfrac(t))

    #---------------------------------------------------------------------------
    # Inner and outer radii as a function of time.
    #---------------------------------------------------------------------------
    def rInner(self, t):
        return self.r0*self.hfrac(t)

    def rOuter(self, t):
        return self.r1*self.hfrac(t)

    #---------------------------------------------------------------------------
    # Inner and outer pressures as a function of time.
    #---------------------------------------------------------------------------
    def Pinner(self, t):
        return self.P0/self.hfrac(t)**(2.0*self.gamma/self.gamma1)

    def Pouter(self, t):
        return self.P1/self.hfrac(t)**(2.0*self.gamma/self.gamma1)

    #---------------------------------------------------------------------------
    # Inner and outer velocities as a function of time.
    #---------------------------------------------------------------------------
    def vrInner(self, t):
        return -self.r0 * t/(self.tau*self.tau * self.hfrac(t))

    def vrOuter(self, t):
        return -self.r1 * t/(self.tau*self.tau * self.hfrac(t))

    #---------------------------------------------------------------------------
    # The initial radius for a given radius and time.
    #---------------------------------------------------------------------------
    def initialRadius(self, t, r):
        ht = self.hfrac(t)
        assert ht > 0.0
        return r/ht

    #---------------------------------------------------------------------------
    # The allowed radial bounds as a function of time.
    #---------------------------------------------------------------------------
    def rRange(self, t):
        hi = self.hfrac(t)
        return self.r0*hi, self.r1*hi

    #---------------------------------------------------------------------------
    # Check that the given r at time t is in the allowed bounds of the shell.
    #---------------------------------------------------------------------------
    def rInBounds(self, t, r):
        rmin, rmax = self.rRange(t)
        return r/rmin >= 0.99999 and r/rmax <= 1.00001

    #---------------------------------------------------------------------------
    # Return r in the expected boundaries.
    #---------------------------------------------------------------------------
    def boundedR(self, t, r):
        rmin, rmax = self.rRange(t)
        if r < rmin:
            return rmin
        elif r > rmax:
            return rmax
        else:
            return r

    #---------------------------------------------------------------------------
    # The initial density and pressure.
    # Note you must pass in the *initial* radius for the point you want!
    #---------------------------------------------------------------------------
    def rhoInitial(self, ri):
        thpt = 1.0/(self.r1**2 - self.r0**2)
        ack = ((self.r1**2 - ri*ri)*thpt*self.rho0**self.gamma1 +
               (ri*ri - self.r0**2)*thpt*self.rho1**self.gamma1)
        acksgn = sgn(ack)
        return acksgn * abs(ack)**(1.0/self.gamma1)

    def Pinitial(self, ri):
        return self.S*(self.rhoInitial(ri))**self.gamma

    #---------------------------------------------------------------------------
    # The density, velocity, pressure, and specific thermal energy profiles as
    # a function of radius and time.
    #---------------------------------------------------------------------------
    def rho(self, t, r):
        t1 = self.tfrac(t)
        hi = self.hfrac(t)
        rho0i = self.rhoInitial(r/hi)
        return rho0i/hi**(2.0/self.gamma1)

    def vr(self, t, r):
        hi = self.hfrac(t)
        return -r*t/(hi*self.tau)**2

    def P(self, t, r):
        hi = self.hfrac(t)
        P0i = self.Pinitial(r/hi)
        return P0i/hi**(2.0*self.gamma/self.gamma1)
        
    def eps(self, t, r):
        return self.P(t, r)/(self.gamma1*self.rho(t, r))
        
    #---------------------------------------------------------------------------
    # The time derivatives of the density, velocity, and pressure.
    #---------------------------------------------------------------------------
    def rhoDot(self, t, r):
        return -2.0*self.hfracDot(t)*self.rho(t, r)/(self.gamma1*self.hfrac(t))

    def vrDot(self, t, r):
        hi = self.hfrac(t)
        return (2.0*r*t*self.hfracDot(t)/(self.tau*self.tau * hi*hi*hi) -
                (r + t*self.vr(t, r))/(self.tau*self.tau * hi*hi))
        #return 2.0* self.tau**2 *r*t*self.hfracDot(t)/self.hfrac(t) - r - t*self.vr(t, r)
        #return -(self.hfracDot(t)*self.vr(t, r) + self.initialRadius(t, r)/(self.tau*self.tau))/self.hfrac(t)

    def Pdot(self, t, r):
        return -2.0*self.gamma*self.hfracDot(t)*self.P(t, r)/(self.gamma1*self.hfrac(t))

    #---------------------------------------------------------------------------
    # The radial gradient of the radial velocity.
    #---------------------------------------------------------------------------
    def DvrDr(self, t, r):
        return -t/(self.tau*self.tau*self.hfrac(t)**2)

    #---------------------------------------------------------------------------
    # The energies in the system as a function of time.
    #---------------------------------------------------------------------------
    def kineticEnergy(self, t):
        return self.K0 * t*t / self.hfrac(t)**((self.gamma + 1.0)/self.gamma1)

    def thermalEnergy(self, t):
        return self.T0 / self.hfrac(t)**((self.gamma + 1.0)/self.gamma1)

    def totalEnergy(self, t):
        return (self.K0 * t*t + self.T0)/self.hfrac(t)**((self.gamma + 1.0)/self.gamma1)

    #---------------------------------------------------------------------------
    # The profiles for the energy computations at time 0.
    #---------------------------------------------------------------------------
    def _K0profile(self, r):
        return 0.5/(self.tau**4)*r*r*self.rhoInitial(r)

    def _T0profile(self, r):
        return self.Pinitial(r)/self.gamma1
