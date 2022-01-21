#===============================================================================
# Class to generate and hold a 2-D node distribution for use in test cases.
#===============================================================================
from math import *
from Spheral import *

import random

class Generate2dTestSetup:

    #===========================================================================
    # Create a set of NodeLists and a DataBase for use in the tests.
    #===========================================================================
    def __init__(self,
                 asph = False,
                 seed = 'random',
                 nodesPerh = 2.01,
                 n1 = 1000,
                 n2 = 2500,
                 n3 = 500,
                 rmin1 = Vector2d(0.0, 0.0),
                 rmax1 = Vector2d(1.0, 1.0),
                 rmin2 = Vector2d(1.0, 0.0),
                 rmax2 = Vector2d(1.5, 1.0),
                 rmin3 = Vector2d(1.5, 0.0),
                 rmax3 = Vector2d(2.0, 1.0)):

        self.n1 = n1
        self.n2 = n2
        self.n3 = n3

        self.g = random.Random()

        self.cache = []

        neighborSearchType = Neighbor2d.NeighborSearchType.GatherScatter
        numGridLevels = 10
        topGridCellSize = 8.0
        origin = Vector2d(0.0, 0.0)
        kernelExtent = 2.0
        hmax = 0.5

        vol1 = (rmax1.x - rmin1.x)*(rmax1.y - rmin1.y)
        vol2 = (rmax2.x - rmin2.x)*(rmax2.y - rmin2.y)
        vol3 = (rmax3.x - rmin3.x)*(rmax3.y - rmin3.y)
        rho1 = 1.0
        rho2 = 1.0
        rho3 = 1.0
        m1 = vol1*rho1 * n1/(n1*n1 + 1e-30)
        m2 = vol2*rho2 * n2/(n2*n2 + 1e-30)
        m3 = vol3*rho3 * n3/(n3*n3 + 1e-30)

        self.eos = GammaLawGasMKS2d(5.0/3.0, 1.0)
        self.WT = TableKernel2d(BSplineKernel2d())

        # Construct the NodeLists to be distributed
        self.dataBase = DataBase2d()
        self.dataBase.updateConnectivityMap()
        if asph:
            self.nodes1 = AsphNodeList2d("nodes 1", self.eos, self.WT, self.WT)
            self.nodes2 = AsphNodeList2d("nodes 2", self.eos, self.WT, self.WT)
            self.nodes3 = AsphNodeList2d("nodes 3", self.eos, self.WT, self.WT)
        else:
            self.nodes1 = SphNodeList2d("nodes 1", self.eos, self.WT, self.WT)
            self.nodes2 = SphNodeList2d("nodes 2", self.eos, self.WT, self.WT)
            self.nodes3 = SphNodeList2d("nodes 3", self.eos, self.WT, self.WT)
        for nodes, n, rmin, rmax, m, rho in ((self.nodes1, n1, rmin1, rmax1, m1, rho1),
                                             (self.nodes2, n2, rmin2, rmax2, m2, rho2),
                                             (self.nodes3, n3, rmin3, rmax3, m3, rho3)):

            if seed == 'random':
                xyNodes = self.randomDistribute(n, rmin, rmax)
            elif seed == 'lattice':
                xyNodes = self.latticeDistribute(n, rmin, rmax)

            nodes.numInternalNodes = n
            nodes.nodesPerSmoothingScale = nodesPerh
            nodes.hmax = hmax
            Hi = self.determineH(n, rmin, rmax, nodesPerh)
            nodes.mass(ScalarField2d("tmp", nodes, m))
            nodes.massDensity(ScalarField2d("tmp", nodes, rho))
            nodes.Hfield(SymTensorField2d("tmp", nodes, Hi))
            nodes.updateWeight(self.dataBase.connectivityMap())
            for i in xrange(n):
                nodes.positions()[i] = xyNodes[i]

            neighbor = NestedGridNeighbor2d(nodes,
                                            neighborSearchType,
                                            numGridLevels,
                                            topGridCellSize,
                                            origin,
                                            kernelExtent)
            nodes.registerNeighbor(neighbor)
            self.cache.append(neighbor)

            self.dataBase.appendNodeList(nodes)

        return
    
    #===========================================================================
    # Calculate one over the smoothing scale for the given number of nodes and
    # volume.
    #===========================================================================
    def determineH(self, nGlobal, rmin, rmax,
                   nNodesPerh = 2.01):
    
        if nGlobal > 0:
            vol = (rmax.y - rmin.y) * (rmax.x  - rmin.x)
            assert vol > 0.0
            dV = vol/nGlobal
            dx = sqrt(dV)
            hi = 1.0/(nNodesPerh*dx)
            Hi = SymTensor2d(hi, 0.0,
                             0.0, hi)
            return Hi

        else:
            return SymTensor2d()

    #===========================================================================
    # Distribute nodes randomly in the given volume.
    #===========================================================================
    def randomDistribute(self,
                         nNodesGlobal,   # global number of nodes in this nodelist
                         rmin, rmax):    # total simulation volume

        nodePositions = []
        for globalNodeID in xrange(nNodesGlobal):
            nodePositions.append(Vector2d(self.g.uniform(rmin.x, rmax.x),
                                          self.g.uniform(rmin.y, rmax.y)))

        assert len(nodePositions) == nNodesGlobal
        return nodePositions

    #===========================================================================
    # Distribute nodes on a lattice in the given volume.
    #===========================================================================
    def latticeDistribute(self,
                          n,              # global number of nodes in this nodelist
                          rmin, rmax):    # total simulation volume

        nodePositions = []

        if n > 0:
            nx = int(sqrt(n) + 1e-5)
            assert nx*nx == n
            dx = (rmax.x - rmin.x)/nx
            dy = (rmax.y - rmin.y)/nx

            for ix in xrange(nx):
                for iy in xrange(nx):
                    nodePositions.append(Vector2d((ix + 0.5)*dx, (iy + 0.5)*dy))

        assert len(nodePositions) == n
        return nodePositions
