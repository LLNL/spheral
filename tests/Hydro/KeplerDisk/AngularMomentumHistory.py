from math import *
import Spheral
import mpi

#-------------------------------------------------------------------------------
# A class for tracking the history of the angular momentum per node.
# Beware : this thing computes and stores a new field for every sample operation.
#          That could get memory intensive eventually!
#-------------------------------------------------------------------------------
class AngularMomentumHistory:

    def __init__(self,
                 nodeList):
        self.restart = Spheral.RestartableObject(self)
        self.nodeList = nodeList
        self.cycleHistory = []
        self.timeHistory = []
        self.LzHistory = []

        # Figure out the dimensionality.
        FieldConstructor = None
        if isinstance(nodeList, Spheral.NodeList1d):
            FieldConstructor = Spheral.ScalarField1d
        elif isinstance(nodeList, Spheral.NodeList2d):
            FieldConstructor = Spheral.ScalarField2d
        elif isinstance(nodeList, Spheral.NodeList3d):
            FieldConstructor = Spheral.ScalarField3d
        assert FieldConstructor is not None

        return

    def sample(self, cycle, t, dt):
        self.cycleHistory.append(cycle)
        self.timeHistory.append(t)

        # Grab the current state.
        pos = self.nodeList.positions()
        vel = self.nodeList.velocity()

        # Build the new Field for the angular momentum.
        self.LzHistory.append(ScalarField("Lz %i %g" % (cycle, t), self.nodeList))
        Lz = self.LzHistory[-1]

        # Walk the nodes and build Lz.
        for i in xrange(self.nodeList.numInternalNodes):
            Lz[i] = pos[i].cross(vel[i]).z

        return

    def label(self):
        return "AngularMomentumHistory"

    def dumpState(self, file, path):
        file.writeObject(self.cycleHistory, path + "/cycleHistory")
        file.writeObject(self.timeHistory, path + "/timeHistory")
        file.writeObject(len(self.LzHistory), path + "/lenLz")
        for i in xrange(len(self.LzHistory)):
            file.write(self.LzHistory[i], path + "/Lz%i" % i)
        return

    def restoreState(self, file, path):
        self.cycleHistory = file.readObject(path + "/cycleHistory")
        self.timeHistory = file.readObject(path + "/timeHistory")
        n = file.readObject(path + "/lenLz")
        for i in xrange(n):
            field = ScalarField("tmp", self.nodeList)
            file.read(field, path + "/Lz%i" % i)
            self.LzHistory.append(field)
        return

