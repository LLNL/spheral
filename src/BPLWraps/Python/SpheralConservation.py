# SpheralConservation
from Spheral import *
import loadmpi
mpi, procID, nProcs = loadmpi.loadmpi()


#-------------------------------------------------------------------------------
# Conservation
#-------------------------------------------------------------------------------
class SpheralConservation(RestartableObject):

    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self, dataBase,
                 packages = []):
        RestartableObject.__init__(self)
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
        self.origin = GeomVector3d(0.0)

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
        for nodeList in self.dataBase.nodeLists():
            for nodeID in xrange(nodeList.numInternalNodes):
                total += nodeList.mass()[nodeID]
        return mpi.allreduce(total, mpi.SUM)

    #---------------------------------------------------------------------------
    # Determine the current total linear momentum.
    #---------------------------------------------------------------------------
    def findTotalPmom(self):
        total = GeomVector3d(0.0)

        # Tally momentum from nodelists.
        for nodeList in self.dataBase.fluidNodeLists():
            for nodeID in xrange(nodeList.numInternalNodes):
                total.x += (nodeList.mass()[nodeID]*
                            nodeList.velocity()[nodeID].x)
                if self.dataBase.nDim > 1:
                    total.y += (nodeList.mass()[nodeID]*
                                nodeList.velocity()[nodeID].y)
                if self.dataBase.nDim > 2:
                    total.z += (nodeList.mass()[nodeID]*
                                nodeList.velocity()[nodeID].z)

        # Tally momentum from packages.
        for package in self.packages:
            packageValue = package.extraMomentum()
            if self.dataBase.nDim == 1:
                total += GeomVector3d(packageValue.x, 0.0, 0.0)
            elif self.dataBase.nDim == 2:
                total += GeomVector3d(packageValue.x, packageValue.y, 0.0)
            else:
                total += packageValue

        return mpi.allreduce(total, mpi.SUM)

    #---------------------------------------------------------------------------
    # Determine the current total angular momentum, with reference to the
    # stored origin.
    #---------------------------------------------------------------------------
    def findTotalAmom(self):
        total = GeomVector3d(0.0)
        for nodeList in self.dataBase.fluidNodeLists():
            for nodeID in xrange(nodeList.numInternalNodes):
                # Find the displacement of this node from the origin.
                dr = GeomVector3d(0.0)
                dr.x = nodeList.positions()[nodeID].x - self.origin.x
                if self.dataBase.nDim > 1:
                    dr.y = nodeList.positions()[nodeID].y - self.origin.y
                if self.dataBase.nDim > 2:
                    dr.z = nodeList.positions()[nodeID].z - self.origin.z

                # Now add this node angular momentum.
                if self.dataBase.nDim == 2:
                    total.z += (nodeList.mass()[nodeID]*
                                (dr.x*nodeList.velocity()[nodeID].y -
                                 dr.y*nodeList.velocity()[nodeID].x))
                elif self.dataBase.nDim == 3:
                    total.x += (nodeList.mass()[nodeID]*
                                (dr.y*nodeList.velocity()[nodeID].z -
                                 dr.z*nodeList.velocity()[nodeID].y))
                    total.y += (nodeList.mass()[nodeID]*
                                (dr.z*nodeList.velocity()[nodeID].x -
                                 dr.x*nodeList.velocity()[nodeID].z))
                    total.z += (nodeList.mass()[nodeID]*
                                (dr.x*nodeList.velocity()[nodeID].y -
                                 dr.y*nodeList.velocity()[nodeID].x))
        return mpi.allreduce(total, mpi.SUM)

    #---------------------------------------------------------------------------
    # Determine the current total kinetic energy.
    #---------------------------------------------------------------------------
    def findTotalKE(self):
        total = 0.0
        for nodeList in self.dataBase.fluidNodeLists():
            for nodeID in xrange(nodeList.numInternalNodes):
                total += nodeList.mass()[nodeID]*nodeList.velocity()[nodeID].magnitude2()
        return 0.5*mpi.allreduce(total, mpi.SUM)

    #---------------------------------------------------------------------------
    # Determine the current total thermal energy.
    #---------------------------------------------------------------------------
    def findTotalTE(self):
        total = 0.0
        for nodeList in self.dataBase.fluidNodeLists():
            for nodeID in xrange(nodeList.numInternalNodes):
                total += nodeList.mass()[nodeID]*nodeList.specificThermalEnergy()[nodeID]
        return mpi.allreduce(total, mpi.SUM)

    #---------------------------------------------------------------------------
    # Determine the current total package (or "external") energy.
    #---------------------------------------------------------------------------
    def findTotalPackageEnergy(self):
        total = 0.0
        for package in self.packages:
            total += package.extraEnergy()
        return mpi.allreduce(total, mpi.SUM)

    #---------------------------------------------------------------------------
    # Write the history to the given file.
    #---------------------------------------------------------------------------
    def writeHistory(self, filename):
        f = open(filename, 'w')
        labels = ['cycle', 'time',
                  'Mass',
                  'Lin Mom Mag', 'Lin Mom X', 'Lin Mom Y', 'Lin Mom Z',
                  'Ang Mom Mag', 'Ang Mom X', 'Ang Mom Y', 'Ang Mom Z',
                  'Total E', 'Kin E', 'Therm E']
        f.write('#')
        for lab in labels:
            f.write('%14s' % lab)
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
                        self.TEHistory[i]]:
                f.write('%14.8g' % var)
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
