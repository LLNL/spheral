from NodeGeneratorBase import *

#-------------------------------------------------------------------------------
# Combine the results of a set of NodeDistribtion generators into one
# distribution.
#-------------------------------------------------------------------------------
class CompositeNodeDistribution(NodeGeneratorBase):

    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self,
                 *generators):
        self._generators = generators

        # Copy the local info from each generator to this one.
        self.positions = []
        self.m = []
        self.rho = []
        self.H = []
        self.globalIDs = []
        offset = 0
        for g in generators:
            for i in range(g.localNumNodes()):
                self.positions.append(g.localPosition(i))
                self.m.append(g.localMass(i))
                self.rho.append(g.localMassDensity(i))
                self.H.append(g.localHtensor(i))
                self.globalIDs.append(offset + g.globalIDs[i])
            offset += g.globalNumNodes()

        NodeGeneratorBase.__init__(self, False, self.positions, self.m, self.rho, self.H)

        return

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        return self.positions[i]

    #---------------------------------------------------------------------------
    # Get the mass for the given node index.
    #---------------------------------------------------------------------------
    def localMass(self, i):
        return self.m[i]

    #---------------------------------------------------------------------------
    # Get the mass density for the given node index.
    #---------------------------------------------------------------------------
    def localMassDensity(self, i):
        return self.rho[i]

    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        return self.H[i]
