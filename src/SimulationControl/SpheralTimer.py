import time
import mpi

################################################################################
# Provide a timing class
################################################################################
from SpheralCompiledPackages import RestartableObject

class SpheralTimer:

    def __init__(self, label=None):
        self.__label = label
        self.numInvocations = 0
        self.lastStartTime = 0.0
        self.lastStopTime = 0.0
        self.lastInterval = 0.0
        self.elapsedTime = 0.0
        self.inTimeCycle = 0
        self.restart = RestartableObject(self)

    def start(self):
        if self.inTimeCycle:
            print("Error in timer, attempt to start before last stop.")
        else:
            self.lastStartTime = time.time()
            self.inTimeCycle = 1

    def stop(self):
        if not self.inTimeCycle:
            print("Error in timer, attempt to stop before start.")
        else:
            self.lastStopTime = time.time()
            self.inTimeCycle = 0
            self.numInvocations = self.numInvocations + 1
            self.lastInterval = self.lastStopTime - self.lastStartTime
            self.elapsedTime = self.elapsedTime + self.lastInterval

    def printStatus(self):
        print("################################################################################")
        if self.__label:
            print("Timing statistics for ", self.__label)
        print("All reports listed as (min, max, avg) across processors.")
        print("Last interval time: \t", self.globalStatistics(self.lastInterval))
        print("Total elapsed time: \t", self.globalStatistics(self.elapsedTime))
        print("Total number of times invoked: \t", self.numInvocations)
        print("Average time interval: \t", self.globalStatistics(self.elapsedTime/(self.numInvocations + 1.0e-30)))
        print("################################################################################")

    def globalStatistics(self, var):
        minVar = mpi.allreduce(var, mpi.MIN)
        maxVar = mpi.allreduce(var, mpi.MAX)
        avgVar = mpi.allreduce(var, mpi.SUM)/mpi.procs
        return minVar, maxVar, avgVar

    def label(self):
        result = "SpheralTimer"
        if self.__label:
            result += "_" + self.__label.replace(" ", "_")
        return result

    def dumpState(self, file, path):
        file.writeObject(self.numInvocations, path + "/numInvocations")
        file.writeObject(self.lastStartTime, path + "/lastStartTime")
        file.writeObject(self.lastStopTime, path + "/lastStopTime")
        file.writeObject(self.lastInterval, path + "/lastInterval")
        file.writeObject(self.elapsedTime, path + "/elapsedTime")
        file.writeObject(self.inTimeCycle, path + "/inTimeCycle")
        file.writeObject(self.__label, path + "/__label")
        return

    def restoreState(self, file, path):
        self.numInvocations = file.readObject(path + "/numInvocations")
        self.lastStartTime = file.readObject(path + "/lastStartTime")
        self.lastStopTime = file.readObject(path + "/lastStopTime")
        self.lastInterval = file.readObject(path + "/lastInterval")
        self.elapsedTime = file.readObject(path + "/elapsedTime")
        self.inTimeCycle = file.readObject(path + "/inTimeCycle")
        self.__label = file.readObject(path + "/__label")
        return

