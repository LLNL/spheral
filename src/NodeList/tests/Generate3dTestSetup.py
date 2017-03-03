#===============================================================================
# Class to generate and hold a 3-D node distribution for use in test cases.
#===============================================================================
from math import *
from Spheral import *

import random

class Generate3dTestSetup:

    #===========================================================================
    # Create a set of NodeLists and a DataBase for use in the tests.
    #===========================================================================
    def __init__(self,
                 asph = False,
                 seed = 'random',
                 nodesPerh = 2.01,
                 nx1 = 10,
                 nx2 = 20,
                 nx3 = 5,
                 rmin1 = Vector3d(0.0, 0.0, 0.0),
                 rmax1 = Vector3d(1.0, 1.0, 1.0),
                 rmin2 = Vector3d(1.0, 0.0, 0.0),
                 rmax2 = Vector3d(1.5, 1.0, 1.0),
                 rmin3 = Vector3d(1.5, 0.0, 0.0),
                 rmax3 = Vector3d(2.0, 1.0, 1.0)):

        n1 = nx1**3
        n2 = nx2**3
        n3 = nx3**3

        self.nx1 = nx1
        self.nx2 = nx2
        self.nx3 = nx3
        self.n1 = n1
        self.n2 = n2
        self.n3 = n3

        self.g = random.Random()

        self.cache = []

        neighborSearchType = Neighbor3d.NeighborSearchType.GatherScatter
        numGridLevels = 10
        topGridCellSize = 8.0
        origin = Vector3d(0.0, 0.0, 0.0)
        kernelExtent = 2.0
        hmax = 0.5

        vol1 = (rmax1.x - rmin1.x)*(rmax1.y - rmin1.y)*(rmax1.z - rmin1.z)
        vol2 = (rmax2.x - rmin2.x)*(rmax2.y - rmin2.y)*(rmax2.z - rmin2.z)
        vol3 = (rmax3.x - rmin3.x)*(rmax3.y - rmin3.y)*(rmax3.z - rmin3.z)
        rho1 = 1.0
        rho2 = 1.0
        rho3 = 1.0
        m1 = vol1*rho1 * n1/(n1*n1 + 1e-30)
        m2 = vol2*rho2 * n2/(n2*n2 + 1e-30)
        m3 = vol3*rho3 * n3/(n3*n3 + 1e-30)

        self.eos = GammaLawGasMKS3d(5.0/3.0, 1.0)

        # Construct the NodeLists to be distributed
        self.dataBase = DataBase3d()
        self.dataBase.updateConnectivityMap()
        if asph:
            self.nodes1 = AsphNodeList3d(self.eos)
            self.nodes2 = AsphNodeList3d(self.eos)
            self.nodes3 = AsphNodeList3d(self.eos)
        else:
            self.nodes1 = SphNodeList3d(self.eos)
            self.nodes2 = SphNodeList3d(self.eos)
            self.nodes3 = SphNodeList3d(self.eos)
        for nodes, nx, rmin, rmax, m, rho in ((self.nodes1, nx1, rmin1, rmax1, m1, rho1),
                                              (self.nodes2, nx2, rmin2, rmax2, m2, rho2),
                                              (self.nodes3, nx3, rmin3, rmax3, m3, rho3)):
            n = nx*nx*nx

            if seed == 'random':
                xyzNodes = self.randomDistribute(n, rmin, rmax)
            elif seed == 'lattice':
                xyzNodes = self.latticeDistribute(nx, rmin, rmax)

            nodes.numInternalNodes = n
            nodes.nodesPerSmoothingScale = nodesPerh
            nodes.hmax = hmax
            Hi = self.determineH(n, rmin, rmax, nodesPerh)
            nodes.setMass(ScalarField3d(nodes, m))
            nodes.setMassDensity(ScalarField3d(nodes, rho))
            nodes.setHfield(SymTensorField3d(nodes, Hi))
            nodes.updateWeight(self.dataBase.connectivityMap())
            for i in xrange(n):
                nodes.positions()[i] = xyzNodes[i]

            neighbor = NestedGridNeighbor3d(nodes,
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
    def determineH(self, n, rmin, rmax,
                   nNodesPerh = 2.01):

        if n > 0:
            vol = (rmax.x  - rmin.x) * (rmax.y - rmin.y) * (rmax.z - rmin.z)
            assert vol > 0.0
            dV = vol/n
            dx = dV**(1.0/3.0)
            hi = 1.0/(nNodesPerh*dx)
            Hi = SymTensor3d(hi, 0.0, 0.0,
                             0.0, hi, 0.0,
                             0.0, 0.0, hi)
            return Hi

        else:
            return SymTensor3d()

    #===========================================================================
    # Distribute nodes randomly in the given volume.
    #===========================================================================
    def randomDistribute(self,
                         n,              # global number of nodes in this nodelist
                         rmin, rmax):    # total simulation volume

        nodePositions = []
        for globalNodeID in xrange(n):
            nodePositions.append(Vector3d(self.g.uniform(rmin.x, rmax.x),
                                          self.g.uniform(rmin.y, rmax.y),
                                          self.g.uniform(rmin.z, rmax.z)))

        assert len(nodePositions) == n
        return nodePositions

    #===========================================================================
    # Distribute nodes on a lattice in the given volume.
    #===========================================================================
    def latticeDistribute(self,
                          nx,             # 
                          rmin, rmax):    # total simulation volume

        nodePositions = []

        n = nx*nx*nx
        if n > 0:
            dx = (rmax.x - rmin.x)/nx
            dy = (rmax.y - rmin.y)/nx
            dz = (rmax.z - rmin.z)/nx

            for ix in xrange(nx):
                for iy in xrange(nx):
                    for iz in xrange(nx):
                        nodePositions.append(Vector3d((ix + 0.5)*dx,
                                                      (iy + 0.5)*dy,
                                                      (iz + 0.5)*dz))

        assert len(nodePositions) == n
        return nodePositions
