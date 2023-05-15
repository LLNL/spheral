#-------------------------------------------------------------------------------
# pnorm
# A class to compute the "p norm" of a given vector of data.  This is defined as
#    ||x||_p = ( \sum_i |x_i|^p )^(1/p)
#
# The special case of p -> \infinity is given by
#    ||x||_\infinity = max( |x_i| )
#
# Optionally we can use the "grid p norm" modification, which takes into account
# resolution as
#    ||x||_gp = ( \sum_i dx (|x_i|)^p )^(1/p)
#-------------------------------------------------------------------------------
import numpy as np
from numpy import linalg as la

class Pnorm:

    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self,
                 vectorData,
                 positionData = None,
                 ansData = None):
        self.positionData = np.array(positionData)
        if ansData is None:
            self.vectorData = np.absolute(np.array(vectorData))
        else:
            self.vectorData = np.absolute(np.array(vectorData) - np.array(ansData))
        return

    #---------------------------------------------------------------------------
    # Compute the slice weighting, i.e., checking if the points are in range.
    #---------------------------------------------------------------------------
    def computeSliceWeighting(self, positionData, rmin, rmax):

        # Make a copy to work on, and sort it.
        n = len(positionData)
        rData = list(zip(positionData[:], list(range(n))))
        rData.sort()

        if rmin is None:
            rmin = min([x[0] for x in rData])
        if rmax is None:
            rmax = max([x[0] for x in rData])

        n = len(positionData)
        weightData = np.zeros(n)
        for i in range(n):
            if positionData[i] >= rmin and positionData[i] <= rmax:
                weightData[i] = 1.0
        return weightData

    #---------------------------------------------------------------------------
    # Compute the grid weighting based on the given position data.
    #---------------------------------------------------------------------------
    def computeGridWeighting(self, positionData, rmin, rmax):

        # Make a copy to work on, and sort it.
        n = len(positionData)
        rData = list(zip(positionData[:], list(range(n))))
        rData.sort()

        if rmin is None:
            rmin = min([x[0] for x in rData])
        if rmax is None:
            rmax = max([x[0] for x in rData])

        # Now build up the grid weighting based on the dr steps.
        weightData = np.zeros(n)
        for j in range(n):
            i = rData[j][1]
            if j == 0:
                r0 = max(rmin, min(rmax, rData[j][0]))
            else:
                r0 = max(rmin, min(rmax, 0.5*(rData[j-1][0] + rData[j][0])))
            if j == n - 1:
                r1 = max(rmin, min(rmax, rData[j][0]))
            else:
                r1 = max(rmin, min(rmax, 0.5*(rData[j][0] + rData[j+1][0])))
            weightData[i] = r1 - r0
            assert weightData[i] >= 0.0

        # That's it, we now have the grid weighting.
        assert len(weightData) == len(positionData)
        assert min(weightData) >= 0.0
        return weightData

    #---------------------------------------------------------------------------
    # Compute the p norm.
    #---------------------------------------------------------------------------
    def pnorm(self, p,
              rmin = None,
              rmax = None):
        weightData = self.computeSliceWeighting(self.positionData, rmin, rmax)

        if p == "inf":
            Ln = la.norm(weightData*self.vectorData, np.inf)
        else:
            Ln = la.norm(weightData*self.vectorData, p)/max(1e-30, sum(weightData))**(1.0/p)
        return Ln

    #---------------------------------------------------------------------------
    # Compute the grid p norm.
    #---------------------------------------------------------------------------
    def gridpnorm(self, p,
                  rmin = None,
                  rmax = None):

        if p == "inf":
            weightData = self.computeSliceWeighting(self.positionData, rmin, rmax)
            Ln = la.norm(weightData*self.vectorData, np.inf)
        else:
            weightData = self.computeGridWeighting(self.positionData, rmin, rmax)
            Ln = la.norm(weightData*self.vectorData, p)/max(1e-30, sum(weightData))**(1.0/p)
        return Ln
