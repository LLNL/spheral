import unittest
from Spheral import *

#-------------------------------------------------------------------------------
# Base class to unit test the ConstantBoundary boundary condition.
#-------------------------------------------------------------------------------
class ConstantBoundaryTest:

    def testApplyBoundary(self):
        assert self.nodes.numInternalNodes == self.n
        assert self.nodes.numGhostNodes == 0
        self.boundary.setGhostNodes(self.nodes)
        self.boundary.applyGhostBoundary(self.nodes.massDensity())
        self.boundary.applyGhostBoundary(self.field)
        assert self.nodes.numGhostNodes == self.nghost
        assert self.boundary.numConstantNodes == self.nghost
        ghostNodes = self.boundary.ghostNodes(self.nodes)
        assert len(ghostNodes) == self.nghost
        for i in ghostNodes:
            r = self.nodes.positions()[i].magnitude()
            assert r > self.rmax and r < self.rbound
            assert abs(self.field[i] + r) < self.tiny
            assert abs(self.nodes.massDensity()[i] - self.rho) < self.tiny

#-------------------------------------------------------------------------------
# 1-D test.
#-------------------------------------------------------------------------------
class ConstantBoundaryTest1d(ConstantBoundaryTest, unittest.TestCase):

    def setUp(self):
        self.tiny = 1.0e-5

        from DistributeNodes import distributeNodes1d
        gamma = 5.0/3.0
        mu = 1.0
        neighborSearchType = Neighbor1d.NeighborSearchType.GatherScatter
        numGridLevels = 10
        topGridCellSize = 0.25
        origin = Vector1d(0.0)
        kernelExtent = 2.0
        self.rho = 1.0
        H1 = SymTensor1d(1.0/0.01)

        self.eos = GammaLawGasMKS1d(gamma, mu)
        self.nodes = SphNodeList1d(self.eos)
        self.neighbor = NestedGridNeighbor1d(self.nodes,
                                             neighborSearchType,
                                             numGridLevels,
                                             topGridCellSize,
                                             origin,
                                             kernelExtent)
        self.nodes.registerNeighbor(self.neighbor)
        self.n = 100
        self.nghost = 20
        self.rmin = 0.0
        self.rmax = 1.0
        self.rbound = 1.2
        distributeNodes1d([(self.nodes, self.n + self.nghost, (self.rmin, self.rbound))])
        self.nodes.setMass(ScalarField1d(self.nodes, 0.5))
        self.nodes.setHfield(SymTensorField1d(self.nodes, H1))
        self.nodes.setMassDensity(ScalarField1d(self.nodes, self.rho))
        constantNodeIDs = vector_of_int()
        for i in xrange(self.n, self.n + self.nghost):
            constantNodeIDs.append(i)
        self.field = ScalarField1d(self.nodes)
        for i in constantNodeIDs:
            self.field[i] = -(self.nodes.positions()[i].magnitude())
        self.boundary = ConstantBoundary1d(self.nodes, constantNodeIDs)
        assert self.boundary.numConstantNodes == self.nghost
        self.nodes.deleteNodes(constantNodeIDs)
        assert self.nodes.numNodes == self.n
        return

#-------------------------------------------------------------------------------
# 2-D test.
#-------------------------------------------------------------------------------
class ConstantBoundaryTest2d(ConstantBoundaryTest, unittest.TestCase):

    def setUp(self):
        self.tiny = 1.0e-5

        from GenerateNodeDistribution2d import GenerateNodeDistribution2d
        from ParMETISDistributeNodes import distributeNodes2d
        gamma = 5.0/3.0
        mu = 1.0
        neighborSearchType = Neighbor2d.NeighborSearchType.GatherScatter
        numGridLevels = 10
        topGridCellSize = 0.25
        origin = Vector2d(0.0)
        kernelExtent = 2.0
        self.rho = 1.0
        seed = "constantDTheta"

        self.eos = GammaLawGasMKS2d(gamma, mu)
        self.nodes = SphNodeList2d(self.eos)
        self.neighbor = NestedGridNeighbor2d(self.nodes,
                                             neighborSearchType,
                                             numGridLevels,
                                             topGridCellSize,
                                             origin,
                                             kernelExtent)
        self.nodes.registerNeighbor(self.neighbor)
        nRadial, nTheta = 50, 50
        nRadialGhost, nThetaGhost = 10, 50
        self.rmin = 0.0
        self.rmax = 1.0
        self.rbound = 1.2
        generator = GenerateNodeDistribution2d(nRadial, nTheta, self.rho, seed,
                                               rmin = self.rmin,
                                               rmax = self.rbound,
                                               nNodePerh = 2.01)
        n1 = generator.globalNumNodes()
        nodeInfo = distributeNodes2d([(self.nodes, n1, generator)])
        self.nodes.setMassDensity(ScalarField2d(self.nodes, self.rho))
        constantNodeIDs = vector_of_int()
        for i in xrange(n1):
            if self.nodes.positions()[i].magnitude() > self.rmax:
                constantNodeIDs.append(i)
        self.nghost = len(constantNodeIDs)
        self.n = self.nodes.numNodes - self.nghost
        self.field = ScalarField2d(self.nodes)
        for i in constantNodeIDs:
            self.field[i] = -(self.nodes.positions()[i].magnitude())
        self.boundary = ConstantBoundary2d(self.nodes, constantNodeIDs)
        assert self.boundary.numConstantNodes == self.nghost
        self.nodes.deleteNodes(constantNodeIDs)
        assert self.nodes.numNodes == self.n
        return

#-------------------------------------------------------------------------------
# 3-D test.
#-------------------------------------------------------------------------------
class ConstantBoundaryTest3d(ConstantBoundaryTest, unittest.TestCase):

    def setUp(self):
        self.tiny = 1.0e-5

        from GenerateNodeDistribution3d import GenerateNodeDistribution3d
        from ParMETISDistributeNodes import distributeNodes3d
        gamma = 5.0/3.0
        mu = 1.0
        neighborSearchType = Neighbor3d.NeighborSearchType.GatherScatter
        numGridLevels = 10
        topGridCellSize = 10.0
        origin = Vector3d(0.0)
        kernelExtent = 2.0
        self.rho = 1.0
        seed = "lattice"

        self.eos = GammaLawGasMKS3d(gamma, mu)
        self.nodes = SphNodeList3d(self.eos)
        self.neighbor = NestedGridNeighbor3d(self.nodes,
                                             neighborSearchType,
                                             numGridLevels,
                                             topGridCellSize,
                                             origin,
                                             kernelExtent)
        self.nodes.registerNeighbor(self.neighbor)
        nx, ny, nz = 20, 20, 20
        nxGhost, nyGhost, nzGhost = 10, 10, 10
        xmin, xmax = (-1.2, -1.2, -1.2), (1.2, 1.2, 1.2)
        self.rmin = 0.0
        self.rmax = 1.0
        self.rbound = 1.2
        generator = GenerateNodeDistribution3d(nx + nxGhost,
                                               ny + nyGhost,
                                               nz + nzGhost, self.rho, seed,
                                               xmin = xmin,
                                               xmax = xmax,
                                               rmin = self.rmin,
                                               rmax = self.rbound,
                                               nNodePerh = 2.01)
        n1 = generator.globalNumNodes()
        nodeInfo = distributeNodes3d([(self.nodes, n1, generator)])
        self.nodes.setMassDensity(ScalarField3d(self.nodes, self.rho))
        constantNodeIDs = vector_of_int()
        for i in xrange(n1):
            if self.nodes.positions()[i].magnitude() > self.rmax:
                constantNodeIDs.append(i)
        self.nghost = len(constantNodeIDs)
        self.n = self.nodes.numNodes - self.nghost
        self.field = ScalarField3d(self.nodes)
        for i in constantNodeIDs:
            self.field[i] = -(self.nodes.positions()[i].magnitude())
        self.boundary = ConstantBoundary3d(self.nodes, constantNodeIDs)
        assert self.boundary.numConstantNodes == self.nghost
        self.nodes.deleteNodes(constantNodeIDs)
        assert self.nodes.numNodes == self.n
        return

if __name__ == "__main__":
    unittest.main()
