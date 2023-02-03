# SpheralConservation
import mpi
from SpheralCompiledPackages import *

#-------------------------------------------------------------------------------
# Conservation
#-------------------------------------------------------------------------------
class SpheralConservation:

    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self, dataBase,
                 packages = []):
        self.restart = RestartableObject(self)
        self.dataBase = dataBase
        self.packages = packages
        self.cycleHistory = []
        self.timeHistory = []
        self.massHistory = []
        self.pmomHistory = []
        self.amomHistory = []
        self.KEHistory = []
        self.TEHistory = []
        self.EEHistory = []
        self.EHistory = []
        self.Vector = eval("Vector%id" % dataBase.nDim)
        self.origin = self.Vector()

        # Start the conservation history
        self.updateHistory()
        return

    #---------------------------------------------------------------------------
    # Add the current state to the history variables.
    #---------------------------------------------------------------------------
    def updateHistory(self, cycle=0, time=0.0):
        self.cycleHistory.append(cycle)
        self.timeHistory.append(time)
        self.massHistory.append(self.findTotalMass())
        self.pmomHistory.append(self.findTotalPmom())
        self.amomHistory.append(self.findTotalAmom())
        self.KEHistory.append(self.findTotalKE())
        self.TEHistory.append(self.findTotalTE())
        self.EEHistory.append(self.findTotalPackageEnergy())
        self.EHistory.append(self.KEHistory[-1] +
                             self.TEHistory[-1] +
                             self.EEHistory[-1])
        return

    #---------------------------------------------------------------------------
    # Determine the current total mass.
    #---------------------------------------------------------------------------
    def findTotalMass(self):
        total = 0.0
        massFL = self.dataBase.globalMass
        for mass in massFL:
            massValues = mass.internalValues()
            total += sum(list(massValues) + [0.0])
        return mpi.allreduce(total, mpi.SUM)

    #---------------------------------------------------------------------------
    # Determine the current total linear momentum.
    #---------------------------------------------------------------------------
    def findTotalPmom(self):
        total = self.Vector()
        massFL = self.dataBase.globalMass
        velocityFL = self.dataBase.globalVelocity
        for (mass, velocity) in zip(massFL, velocityFL):
            massValues = mass.internalValues()
            velocityValues = velocity.internalValues()
            for mi, vi in zip(massValues, velocityValues):
                total += mi*vi

        # Tally momentum from packages.
        for package in self.packages:
            packageValue = package.extraMomentum()
            total += packageValue

        return mpi.allreduce(total, mpi.SUM)

    #---------------------------------------------------------------------------
    # Determine the current total angular momentum, with reference to the
    # stored origin.
    #---------------------------------------------------------------------------
    def findTotalAmom(self):
        total = Vector3d()
        massFL = self.dataBase.globalMass
        positionFL = self.dataBase.globalPosition
        velocityFL = self.dataBase.globalVelocity

        for (mass, position, velocity) in zip(massFL, positionFL, velocityFL):
            massValues = mass.internalValues()
            positionValues = position.internalValues()
            velocityValues = velocity.internalValues()
            for (mi, ri, vi) in zip(massValues, positionValues, velocityValues):

                # Find the displacement of this node from the origin.
                dr = ri - self.origin

                # Now add this node angular momentum.
                if self.dataBase.nDim == 2:
                    total.z += mi*(dr.x*vi.y - dr.y*vi.x)
                elif self.dataBase.nDim == 3:
                    total += mi * dr.cross(vi)

        return mpi.allreduce(total, mpi.SUM)

    #---------------------------------------------------------------------------
    # Determine the current total kinetic energy.
    #---------------------------------------------------------------------------
    def findTotalKE(self):
        total = 0.0
        massFL = self.dataBase.globalMass
        velocityFL = self.dataBase.globalVelocity

        for (mass, velocity) in zip(massFL, velocityFL):
            massValues = mass.internalValues()
            velocityValues = velocity.internalValues()
            total += sum([mi*vi.magnitude2() for (mi, vi) in zip(massValues, velocityValues)] + [0.0])

        return 0.5*mpi.allreduce(total, mpi.SUM)

    #---------------------------------------------------------------------------
    # Determine the current total thermal energy.
    #---------------------------------------------------------------------------
    def findTotalTE(self):
        total = 0.0
        massFL = self.dataBase.fluidMass
        epsFL = self.dataBase.fluidSpecificThermalEnergy

        for (mass, eps) in zip(massFL, epsFL):
            massValues = mass.internalValues()
            epsValues = eps.internalValues()
            total += sum([mi*epsi for (mi, epsi) in zip(list(mass.internalValues()),
                                                        list(eps.internalValues()))] + [0.0])

        return mpi.allreduce(total, mpi.SUM)

    #---------------------------------------------------------------------------
    # Determine the current total package (or "external") energy.
    #---------------------------------------------------------------------------
    def findTotalPackageEnergy(self):
        total = 0.0
        for package in self.packages:
            total += package.extraEnergy()
        return total  # Note we assume this has already been parallel summed.

    #---------------------------------------------------------------------------
    # Write the history to the given file.
    #---------------------------------------------------------------------------
    def writeHistory(self, filename):
        f = open(filename, 'w')
        labels = ['"cycle"', '"time"',
                  '"Mass"',
                  '"Lin Mom Mag"', '"Lin Mom X"', '"Lin Mom Y"', '"Lin Mom Z"',
                  '"Ang Mom Mag"', '"Ang Mom X"', '"Ang Mom Y"', '"Ang Mom Z"',
                  '"Total E"', '"Kin E"', '"Therm E"', '"Pkg E"']
        f.write('#')
        for lab in labels:
            f.write('%14s ' % lab)
        f.write('\n')
        for i in xrange(len(self.cycleHistory)):
            for var in [self.cycleHistory[i], self.timeHistory[i],
                        self.massHistory[i],
                        self.pmomHistory[i].magnitude(),
                        self.pmomHistory[i].x,
                        self.pmomHistory[i].y,
                        self.pmomHistory[i].z,
                        self.amomHistory[i].magnitude(),
                        self.amomHistory[i].x,
                        self.amomHistory[i].y,
                        self.amomHistory[i].z,
                        self.EHistory[i],
                        self.KEHistory[i],
                        self.TEHistory[i],
                        self.EEHistory[i]]:
                f.write('%14.8g ' % var)
            f.write('\n')
        f.close()
        return

    #---------------------------------------------------------------------------
    # label
    #---------------------------------------------------------------------------
    def label(self):
        return "SpheralConservation"

    #---------------------------------------------------------------------------
    # dumpState
    #---------------------------------------------------------------------------
    def dumpState(self, file, path):
        file.writeObject(self.cycleHistory, path + "/cycleHistory")
        file.writeObject(self.timeHistory, path + "/timeHistory")
        file.writeObject(self.massHistory, path + "/massHistory")
        file.writeObject(self.pmomHistory, path + "/pmomHistory")
        file.writeObject(self.amomHistory, path + "/amomHistory")
        file.writeObject(self.KEHistory, path + "/KEHistory")
        file.writeObject(self.TEHistory, path + "/TEHistory")
        file.writeObject(self.EEHistory, path + "/EEHistory")
        file.writeObject(self.EHistory, path + "/EHistory")
        file.writeObject(self.origin, path + "/origin")

    #---------------------------------------------------------------------------
    # restoreState
    #---------------------------------------------------------------------------
    def restoreState(self, file, path):
        self.cycleHistory = file.readObject(path + "/cycleHistory")
        self.timeHistory = file.readObject(path + "/timeHistory")
        self.massHistory = file.readObject(path + "/massHistory")
        self.pmomHistory = file.readObject(path + "/pmomHistory")
        self.amomHistory = file.readObject(path + "/amomHistory")
        self.KEHistory = file.readObject(path + "/KEHistory")
        self.TEHistory = file.readObject(path + "/TEHistory")
        self.EEHistory = file.readObject(path + "/EEHistory")
        self.EHistory = file.readObject(path + "/EHistory")
        self.origin = file.readObject(path + "/origin")
