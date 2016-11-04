from math import *
import Spheral
import mpi

#-------------------------------------------------------------------------------
# A class for tracking the history of a given set of nodes.
#-------------------------------------------------------------------------------
class NodeHistory:

    def __init__(self,
                 nodeList,
                 nodeIndicies,
                 sampleMethod,
                 filename,
                 header = None,
                 labels = None):
        self.restart = Spheral.RestartableObject(self)
        self.nodeList = nodeList
        self.nodeIndicies = nodeIndicies
        self.sampleMethod = sampleMethod
        self.filename = filename
        self.cycleHistory = []
        self.timeHistory = []
        self.sampleHistory = []

        # Figure out the dimensionality.
        FieldConstructor = None
        if isinstance(nodeList, Spheral.NodeList1d):
            FieldConstructor = Spheral.IntField1d
        elif isinstance(nodeList, Spheral.NodeList2d):
            FieldConstructor = Spheral.IntField2d
        elif isinstance(nodeList, Spheral.NodeList3d):
            FieldConstructor = Spheral.IntField3d
        assert FieldConstructor is not None

        # Store the set of nodes we're going to sample as a field of flags.
        # This should automatically be safe as NodeLists/Fields get renumbered,
        # redistributed, deleted, added, or what have you.
        self.nodeFlags = FieldConstructor("flag nodes", nodeList, 0)
        if type(self.nodeIndicies) == list:
            for i in nodeIndicies:
                assert i >= 0 and i < nodeList.numInternalNodes
                self.nodeFlags[i] = 1
        else:
            self.currentNodeIndicies()

        # Open the history file.
        self.file = None
        if mpi.rank == 0:
            self.file = open(self.filename, "w")
            assert self.file is not None

            # Write the optional header string.
            if header:
                self.file.write(header + "\n")

            # Write the optional label line
            if labels:
                self.file.write(("# " + ((len(labels) + 2)*'"%20s" ') + "\n") % (("cycle", "time") + labels))

        return

    def currentNodeIndicies(self):
        if type(self.nodeIndicies) == list:
            return [i for i in range(self.nodeList.numInternalNodes)
                    if self.nodeFlags[i] == 1]
        else:
            result = self.nodeIndicies(self.nodeList)
            self.nodeFlags.Zero()
            for i in result:
                assert i >= 0 and i < self.nodeList.numInternalNodes
                self.nodeFlags[i] = 1
            return result

    def sample(self, cycle, t, dt):

        # Get the set of nodes.
        nodeIndicies = self.currentNodeIndicies()

        # Get the result of the sampling method.
        result = self.sampleMethod(self.nodeList, nodeIndicies)

        # Update our history variables.
        self.cycleHistory.append(cycle)
        self.timeHistory.append(t)
        self.sampleHistory.append(result)

        # Update the history file.
        if mpi.rank == 0:
            assert self.file is not None
            if isinstance(result, tuple):
                samplestr = ""
                for x in result:
                    samplestr += str(x) + " "
            else:
                samplestr = str(result)
            self.file.write("%i \t %g \t %s\n" % (cycle, t, samplestr))
            self.file.flush()

        return

    def flushHistory(self):
        if mpi.rank == 0:
            assert self.file is not None
            n = len(self.cycleHistory)
            assert len(self.timeHistory) == n
            assert len(self.sampleHistory) == n
            if mpi.rank == 0:
                for i in xrange(n):
                    if isinstance(self.sampleHistory[i], tuple):
                        samplestr = ""
                        for x in self.sampleHistory[i]:
                            samplestr += str(x) + " "
                    else:
                        samplestr = str(self.sampleHistory[i])
                    self.file.write("%i \t %g \t %s\n" % (self.cycleHistory[i],
                                                          self.timeHistory[i],
                                                          samplestr))
            self.file.flush()
        return
                
    def label(self):
        return "NodeHistory"

    def dumpState(self, file, path):
        file.writeObject(self.filename, path + "/filename")
        file.writeObject(self.cycleHistory, path + "/cycleHistory")
        file.writeObject(self.timeHistory, path + "/timeHistory")
        file.writeObject(self.sampleHistory, path + "/sampleHistory")
        file.write(self.nodeFlags, path + "/nodeFlags")
        return

    def restoreState(self, file, path):
        self.filename = file.readObject(path + "/filename")
        self.cycleHistory = file.readObject(path + "/cycleHistory")
        self.timeHistory = file.readObject(path + "/timeHistory")
        self.sampleHistory = file.readObject(path + "/sampleHistory")
        file.read(self.nodeFlags, path + "/nodeFlags")
        return

