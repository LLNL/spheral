# Find the mix length for the two fluids in a Rayleigh-Taylor calculation.
# This class is written as functor appropriate to be called by the NodeHistory Spheral worker class
# for periodically collecting node time history diagnostics in Spheral run.

import mpi
import numpy as np
import Spheral2d

class RTMixLength:

    def __init__(self, nodes, y0, 
                 percentile = 95.0):  # percent of material above/below where we're measuring the mix length.
        assert percentile >= 0.0 and percentile <= 100.0
        self.restart = Spheral2d.RestartableObject(self)
        self.nodes = nodes
        self.y0 = y0
        self.percentile = percentile

        # We build nodeFlags so that:
        #   pos.y < y0 ==> 0
        #   pos.y > y0 ==> 1
        self.nodeFlags = Spheral2d.IntField("original node distribution", nodes, 0)
        pos = nodes.positions()
        for i in xrange(nodes.numInternalNodes):
            if pos[i].y > y0:
                self.nodeFlags[i] = 1
        return

    def __call__(self, nodes, nodeIndices):
        # Note we're going to ignore the nodeIndices handed in here.
        pos = nodes.positions()
        y0 = np.array(mpi.reduce([pos[i].y for i in xrange(nodes.numInternalNodes) if self.nodeFlags[i] == 0], mpi.SUM, 0))
        y1 = np.array(mpi.reduce([pos[i].y for i in xrange(nodes.numInternalNodes) if self.nodeFlags[i] == 1], mpi.SUM, 0))
        ylow, yhigh = 0.0, 0.0
        if mpi.rank == 0:
            yhigh = np.percentile(y0, self.percentile)
            ylow = np.percentile(y1, 100.0 - self.percentile)
        ylow = mpi.bcast(ylow, 0)
        yhigh = mpi.bcast(yhigh, 0)
        return yhigh, ylow, yhigh - ylow

    def label(self):
        return "RTMixLength"

    def dumpState(self, file, path):
        file.write(self.nodeFlags, path + "/nodeFlags")
        return

    def restoreState(self, file, path):
        file.read(self.nodeFlags, path + "/nodeFlags")
        return

