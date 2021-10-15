from math import *
import mpi

#-------------------------------------------------------------------------------
# Create a physics package to dump excessive velocities into thermal energy.
#-------------------------------------------------------------------------------
class VelocityKiller:

    def __init__(self, vmax, restoreEnergy, dataBase):
        self.vmax = vmax
        self.restoreEnergy = restoreEnergy
        self.vmax2 = vmax*vmax
        self.dataBase = dataBase
        return

    def doit(self, cycle, atime, dt):
        velocities = self.dataBase.fluidVelocity
        specificEnergies = self.dataBase.fluidSpecificThermalEnergy
        assert len(velocities) == len(specificEnergies)
        nthump = 0
        for vel, eps in zip(velocities, specificEnergies):
            assert vel.numInternalElements == eps.numInternalElements
            for i in xrange(vel.numInternalElements):
                vmag2 = vel[i].magnitude2()
                if vmag2 > self.vmax2:
                    deps = 0.5*(vmag2 - self.vmax2)
                    assert deps > 0.0
                    vel[i] *= self.vmax/sqrt(vmag2)
                    if self.restoreEnergy:
                        eps[i] += deps
                    nthump += 1
        print "Velocity killer subdued %i nodes." % mpi.allreduce(nthump, mpi.SUM)
        return

