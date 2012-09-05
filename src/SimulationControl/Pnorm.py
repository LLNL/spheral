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

class Pnorm:

    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self,
                 vectorData,
                 positionData = None):
        self.vectorData = vectorData
        self.positionData = positionData
        return

    #---------------------------------------------------------------------------
    # Compute the grid weighting based on the given position data.
    #---------------------------------------------------------------------------
    def computeGridWeighting(self, positionData, rmin, rmax):

        # Make a copy to work on, and sort it.
        n = len(positionData)
        rData = zip(positionData[:], range(n))
        rData.sort()

        if rmin is None:
            rmin = min(rData)
        if rmax is None:
            rmax = max(rData)

        # Now build up the grid weighting based on the dr steps.
        weightData = [0.0]*n
        for j in xrange(n):
            i = rData[j][1]
            jm = max(j - 1, 0)
            jp = min(j + 1, n - 1)
            rm = max(rmin, min(rmax, rData[jm][0]))
            rp = max(rmin, min(rmax, rData[jp][0]))
            weightData[i] = 0.5*(rp - rm)
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
        weightData = [1.0]*len(self.vectorData)
        Ln, nused = self._pnorm(p,
                                self.vectorData,
                                weightData,
                                rmin,
                                rmax)
        if p == "inf":
            return Ln
        else:
            return Ln*float(nused)**(1.0/p)

    #---------------------------------------------------------------------------
    # Compute the grid p norm.
    #---------------------------------------------------------------------------
    def gridpnorm(self, p,
                  rmin = None,
                  rmax = None):
        weightData = self.computeGridWeighting(self.positionData, rmin, rmax)
        Ln, nused = self._pnorm(p,
                                self.vectorData,
                                weightData,
                                rmin,
                                rmax)
        return Ln

    #---------------------------------------------------------------------------
    # Compute the average p norm.
    #---------------------------------------------------------------------------
    def pnormAverage(self, p,
                     rmin = None,
                     rmax = None):
        weightData = [1.0]*len(self.vectorData)
        Ln, nused = self._pnorm(p,
                                self.vectorData,
                                weightData,
                                rmin,
                                rmax)
        return Ln

    #---------------------------------------------------------------------------
    # Internal method to do the actual computation.
    #---------------------------------------------------------------------------
    def _pnorm(self, p, vectorData, weightData, rmin, rmax):

        # Pre-conditions.
        assert len(vectorData) == len(weightData)
        assert p != 0.0
        assert ((rmin is None) and (rmax is None)) or (len(self.positionData) == len(vectorData))

        # Compute that norm.
        nused = 0
        totalWeight = 0.0
        result = 0.0
        for i in xrange(len(vectorData)):
            if (((rmin is None) or (self.positionData[i] >= rmin)) and
                ((rmax is None) or (self.positionData[i] <= rmax))):
                vd = vectorData[i]
                wd = weightData[i]
                nused += 1
                totalWeight += wd
                if p == "inf":
                    result = max(result, abs(vd))
                else:
                    result += wd*(abs(vd))**p

        if p != "inf":
            assert totalWeight > 0.0
            result = (result/totalWeight)**(1.0/p)

        return result, nused
