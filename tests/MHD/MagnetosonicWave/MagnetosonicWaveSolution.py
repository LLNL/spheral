from math import *

#-------------------------------------------------------------------------------
# The analytic answer for the magnetosonic wave problem.
#-------------------------------------------------------------------------------
class MagnetosonicWaveSolution:

    # Constructor.
    def __init__(self, eos, rho0, rho1, u0, B0, k, h0, waveType = 'fast'):
        self.eos = eos
        self.cs = cs
        self.rho0 = rho0
        self.rho1 = rho1
        self.u0 = u0
        self.B0 = rho
        self.k = k
        self.h0 = h0
        self.waveType = waveType
        return

    # Compute the sound speed.
    def soundSpeed(self):
        Cs = eos.soundSpeed(self.rho0, self.u0)
        B0 = self.B0
        B02 = B0.dot(B0)
        rho0 = self.rho0
        if options.problem == 'fast':
            return sqrt(0.5*(Cs*Cs + B02/rho0 + \
                        sqrt((Cs*Cs + B02/rho0)**2 - \
                              4.0*Cs*Cs*B0.x*B0.x/rho0)))
        else:
            return sqrt(0.5*(Cs*Cs + B02/rho0 - \
                        sqrt((Cs*Cs + B02/rho0)**2 - \
                             4.0*Cs*Cs*B0.x*B0.x/rho0)))

    # Compute and return the solution on the given positions.
    def solution(self, time, xvals):
        omegat = self.k*self.cs*time
        cs = self.soundSpeed()
        rho = [self.rho0*(1.0 + self.rho1*sin(self.k*x - omegat)) for x in xvals]

        # Figure out the velocity perterbation.
        rho0 = self.rho0
        vx = cs
        Bx, By, Bz = options.B0.x, options.B0.y, options.B0.z
        W = cs*cs - Bx*Bx/rho0
        vy = -Bx*By*vx/(rho*W)
        vz = -Bx*Bz*vx/(rho*W)
        v = [Vector3d(vx * sin(self.k*x.x - omegat), \
                      vy * sin(self.k*x.x - omegat), \
                      vz * sin(self.k*x.x - omegat)) for x in xvals]

        # Figure out the magnetic perterbation.
        by = cs*By*vx/W
        bz = cs*Bz*vx/W
        B = [Vector3d(Bx, \
                      by * sin(self.k*x.x - omegat), \
                      bz * sin(self.k*x.x - omegat)) for x in xvals]

        # Figure out the pressure and specific thermal energy.
        P = [self.eos.pressure(rhoi, ui) for (rhoi, ui) in zip(rho, u)]
        if isinstance(self.eos, GammaLawGasMKS3d):
           u = [P[i]/(eos.getGamma()-1)*rho[i] for i in xrange(len(xvals))]
        else:
           u = [self.u0]*len(xvals)

        # And, of course, the SPH smoothing scale.
        h = [self.h0*rhoi/self.rho0 for rhoi in rho]

        return xvals, rho, v, u, B, P, h
