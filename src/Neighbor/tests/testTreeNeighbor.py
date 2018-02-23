#ATS:test(SELF, np=1, label="TreeNeighbor unit tests")
from math import *
import unittest
import random

from Spheral import *

#===============================================================================
# Apply a random rotation operation to the given symmetric tensor.
#===============================================================================
def applyRandomRotation(t, nDimensions):
    if nDimensions == 2:
        theta = random.uniform(0.0, pi)
        R = Tensor2d(cos(theta), -sin(theta),
                     sin(theta), cos(theta))
        t.rotationalTransform(R)
    elif nDimensions == 3:
        theta1 = random.uniform(0.0, pi)
        theta2 = random.uniform(0.0, pi)
        theta3 = random.uniform(0.0, pi)
        R1 = Tensor3d(cos(theta1), -sin(theta1),  0.0,
                      sin(theta1),  cos(theta1),  0.0,
                      0.0,          0.0,          1.0)
        R2 = Tensor3d(cos(theta2),  0.0,        -sin(theta2),
                      0.0,          1.0,          0.0,
                      sin(theta2),  0.0,          cos(theta2))
        R3 = Tensor3d(1.0,          0.0,          0.0,
                      0.0,          cos(theta3), -sin(theta3),
                      0.0,          sin(theta3),  cos(theta3))
        R = R1*R2*R3
        t.rotationalTransform(R)
    return

#===============================================================================
# Base class for the tests, providing the generic set up methods.
#===============================================================================
class SetupNodeDistributions:

    #---------------------------------------------------------------------------
    # Distribute nodes randomly amongst domains.
    #---------------------------------------------------------------------------
    def randomDistribute(self,
                         n,                  # global number of nodes
                         range,              # total simulation volume
                         Vector,             # the vector type
                         SymTensor,          # tensor type
                         nNodesPerh = 2.01): # duh!
        assert Vector.nDimensions == SymTensor.nDimensions
        assert len(range) == Vector.nDimensions
        assert min([len(range[i]) == 2 for i in xrange(Vector.nDimensions)])

        adim = 1.0/SymTensor.nDimensions
        vol = 1.0
        for i in xrange(SymTensor.nDimensions):
            vol *= range[i][1] - range[i][0]
        assert vol > 0.0
        dx0 = (vol/n)**adim

        globalNodeIDs = []
        nodePositions = []
        H = []
        for i in xrange(n):
            globalNodeIDs.append(i)
            args = (random.uniform(range[j][0], range[j][1])
                    for j in xrange(Vector.nDimensions))
            nodePositions.append(Vector(*args))

            dx = [nNodesPerh * random.uniform(0.5, 2.0)*dx0
                  for j in xrange(SymTensor.nDimensions)]
            assert len(dx) == SymTensor.nDimensions
            assert min(dx) > 0.0
            args = [0.0] * SymTensor.nDimensions**2
            for j in xrange(SymTensor.nDimensions):
                args[j + j*SymTensor.nDimensions] = 1.0/dx[j]
            Hi = SymTensor(*tuple(args))
            applyRandomRotation(Hi, SymTensor.nDimensions)
            H.append(Hi)

        assert len(globalNodeIDs) == n
        assert len(nodePositions) == n
        assert len(H) == n
        return globalNodeIDs, nodePositions, H

    #---------------------------------------------------------------------------
    # setUp routine common to all instances.
    #---------------------------------------------------------------------------
    def genericSetUp(self, n1, n2, n3,
                     range1, range2, range3,
                     EOS,
                     makeFluidNodeList,
                     TableKernel,
                     BSplineKernel,
                     Vector,
                     SymTensor,
                     TreeNeighbor,
                     DataBase,
                     Seeder):
        
        # Construct the NodeLists to be distributed
        self.eos = EOS(2.0, 2.0)
        self.WT = TableKernel(BSplineKernel(), 100)
        self.nodes1 = makeFluidNodeList("nodes 1", self.eos, NeighborType=TreeNeighbor)
        self.nodes2 = makeFluidNodeList("nodes 2", self.eos, NeighborType=TreeNeighbor)
        self.nodes3 = makeFluidNodeList("nodes 3", self.eos, NeighborType=TreeNeighbor)
        for nodes, nGlobal, range in ((self.nodes1, n1, range1),
                                      (self.nodes2, n2, range2),
                                      (self.nodes3, n3, range3)):
            globalIDs, xyNodes, H = Seeder(nGlobal,
                                           range,
                                           Vector,
                                           SymTensor)
            n = len(globalIDs)
            nodes.numInternalNodes = n
            for i in xrange(n):
                nodes.mass()[i] = 1.0
                nodes.positions()[i] = xyNodes[i]
                nodes.Hfield()[i] = H[i]

        # Put the distributed NodeLists into a DataBase.
        self.dataBase = DataBase()
        self.dataBase.appendNodeList(self.nodes1)
        self.dataBase.appendNodeList(self.nodes2)
        self.dataBase.appendNodeList(self.nodes3)

        # Update all the Neighbor info.
        #TreeNeighbor.setBoundingBox()
        for n in (self.nodes1, self.nodes2, self.nodes3):
            n.neighbor().updateNodes()

        return

    #---------------------------------------------------------------------------
    # Method called after test is completed.
    #---------------------------------------------------------------------------
    def tearDown(self):
        return

#===============================================================================
# Base class to implement TreeNeighbor test.
#===============================================================================
class TestTreeNeighborBase(SetupNodeDistributions):

    #---------------------------------------------------------------------------
    # The actual test itself!
    #---------------------------------------------------------------------------
    def testThemNeighbors(self):
        from SpheralTestUtilities import findNeighborNodes, checkNeighbors
        import time

        # Iterate over the NodeLists.
        for nodes in self.dataBase.nodeLists():
            pos = nodes.positions()
            H = nodes.Hfield()

            # Randomly select nodes from each NodeList to explicitly test.
            for nodeID in random.sample(range(nodes.numInternalNodes - 1), self.ncheck):
                ri = pos[nodeID]
                Hi = H[nodeID]

                # Have the neighbor objects select neighbors for this node.
                t0 = time.time()
                masterLists = vector_of_vector_of_int()
                coarseNeighbors = vector_of_vector_of_int()
                refineNeighbors = vector_of_vector_of_int()
                self.dataBase.setMasterNodeLists(ri, Hi, masterLists, coarseNeighbors)
                self.dataBase.setRefineNodeLists(ri, Hi, coarseNeighbors, refineNeighbors)
                neighborIDs = []
                offset = 0
                for inds, nds in enumerate(self.dataBase.nodeLists()):
                    neighborIDs.extend([i + offset for i in refineNeighbors[inds]])
                    offset += nds.numInternalNodes
                t1 = time.time()

                # Now build the checks.
                answerIDs = []
                offset = 0
                for nds in self.dataBase.nodeLists():
                    answerIDs.extend([i + offset for i in findNeighborNodes(ri, self.kernelExtent, nds)])
                    offset += nds.numInternalNodes
                t2 = time.time()

                # Check the answer.
                test = checkNeighbors(neighborIDs, answerIDs)
                t3 = time.time()
                if not test:
                    neighborIDs.sort()
                    answerIDs.sort()
                    print 'SPH Tree Neighbor test FAILED'
                    print ' refine: ', neighborIDs
                    print ' answer: ', answerIDs
                    missing = [i for i in answerIDs if i not in neighborIDs]
                    print 'missing: ', missing
                    print 'deltas: ', [((Hi*(pos[i] - ri)).x, (H[i]*(pos[i] - ri)).x) for i in missing]
                else:
                    print "Passed for node %i : %f %f %f" % (nodeID, t1 - t0, t2 - t1, t3 - t2)
                assert test

#===============================================================================
# Radom node distribution -- 1-D.
#===============================================================================
class TestTreeNeighborRandom1d(unittest.TestCase, TestTreeNeighborBase):

    #---------------------------------------------------------------------------
    # Set up method called before test is run.
    #---------------------------------------------------------------------------
    def setUp(self):

        print "--------------------------------------------------------------------------------"
        print "1-D TreeNeighbor random test."
        print "--------------------------------------------------------------------------------"

        self.ncheck = 10

        # Generic parameters for 1-D tests.
        n1 = 1000
        n2 = 2500
        n3 = 500

        range1 = ((-2.0, -1.0),)
        range2 = ((-1.0, 0.5),)
        range3 = ((0.5, 2.0),)

        searchType = GatherScatter
        numGridLevels = 20
        topGridCellSize = 100.0
        self.kernelExtent = 2.0

        self.genericSetUp(n1, n2, n3,
                          range1, range2, range3,
                          GammaLawGasMKS1d,
                          makeFluidNodeList1d,
                          TableKernel1d,
                          BSplineKernel1d,
                          Vector1d,
                          SymTensor1d,
                          TreeNeighbor1d,
                          DataBase1d,
                          self.randomDistribute)

        return

#===============================================================================
# Radom node distribution -- 2-D.
#===============================================================================
class TestTreeNeighborRandom2d(unittest.TestCase, TestTreeNeighborBase):

    #---------------------------------------------------------------------------
    # Set up method called before test is run.
    #---------------------------------------------------------------------------
    def setUp(self):

        print "--------------------------------------------------------------------------------"
        print "2-D TreeNeighbor random test."
        print "--------------------------------------------------------------------------------"

        self.ncheck = 10

        # Generic parameters for 2-D tests.
        n1 = 10000
        n2 = 25000
        n3 = 5000

        range1 = ((-2.0, -1.0), (0.0, 1.0))
        range2 = ((-1.0, -0.5), (0.0, 1.0))
        range3 = ((-0.5, 0.0), (0.0, 1.0))

        self.kernelExtent = 2.0

        self.genericSetUp(n1, n2, n3,
                          range1, range2, range3,
                          GammaLawGasMKS2d,
                          makeFluidNodeList2d,
                          TableKernel2d,
                          BSplineKernel2d,
                          Vector2d,
                          SymTensor2d,
                          TreeNeighbor2d,
                          DataBase2d,
                          self.randomDistribute)

        return

#===============================================================================
# Radom node distribution -- 3-D.
#===============================================================================
class TestTreeNeighborRandom3d(unittest.TestCase, TestTreeNeighborBase):

    #---------------------------------------------------------------------------
    # Set up method called before test is run.
    #---------------------------------------------------------------------------
    def setUp(self):

        print "--------------------------------------------------------------------------------"
        print "3-D TreeNeighbor random test."
        print "--------------------------------------------------------------------------------"

        self.ncheck = 10

        # Generic parameters for 3-D tests.
        n1 = 1000
        n2 = 2500
        n3 = 1500

        range1 = ((0.0, 1.0), (0.0, 1.0), (0.0, 1.0))
        range2 = ((1.0, 1.5), (0.0, 1.0), (0.0, 1.0))
        range3 = ((1.5, 2.0), (0.0, 1.0), (0.0, 1.0))

        searchType = GatherScatter
        numGridLevels = 20
        topGridCellSize = 100.0
        origin = Vector3d(0.0, 0.0, 0.0)
        self.kernelExtent = 2.0

        self.genericSetUp(n1, n2, n3,
                          range1, range2, range3,
                          GammaLawGasMKS3d,
                          makeFluidNodeList3d,
                          TableKernel3d,
                          BSplineKernel3d,
                          Vector3d,
                          SymTensor3d,
                          TreeNeighbor3d,
                          DataBase3d,
                          self.randomDistribute)


        return

#===============================================================================
# Cylindrical node distribution -- 2-D.
#===============================================================================
class TestTreeNeighborCylindrical2d(unittest.TestCase, TestTreeNeighborBase):

    #---------------------------------------------------------------------------
    # Set up method called before test is run.
    #---------------------------------------------------------------------------
    def setUp(self):

        print "--------------------------------------------------------------------------------"
        print "2-D TreeNeighbor regular cylindrical test."
        print "--------------------------------------------------------------------------------"

        self.ncheck = 50

        from GenerateNodeDistribution2d import GenerateNodeDistribution2d
        from DistributeNodes import distributeNodes2d
        self.eos = GammaLawGasMKS2d(2.0, 2.0)
        self.WT = TableKernel2d(BSplineKernel2d(), 100)
        self.nodes1 = makeFluidNodeList2d("cylindrical nodes 1", self.eos, NeighborType=TreeNeighbor2d)
        self.kernelExtent = 2.0
        generator = GenerateNodeDistribution2d(nRadial = 100,
                                               nTheta = 100,
                                               rho = 1.0,
                                               distributionType = "constantDTheta",
                                               rmin = 0.0,
                                               rmax = 1.0,
                                               theta = 0.5*pi,
                                               nNodePerh = 2.01)
        distributeNodes2d((self.nodes1, generator))
        self.dataBase = DataBase2d()
        self.dataBase.appendNodeList(self.nodes1)

        return

#===============================================================================
# Run the tests
#===============================================================================
if __name__ == "__main__":
    unittest.main()
