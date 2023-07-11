import SolidSpheral
import random
from math import *
import mpi

#-------------------------------------------------------------------------------
# Define a method for seeding flaws.
#-------------------------------------------------------------------------------
def flawDistribution(seed,
                     kWeibull,
                     mWeibull,
                     volume,
                     nodeList,
                     minFlawsPerNode,
                     minTotalFlaws,
                     VectorDoubleField,
                     numGlobalNodes,
                     globalNodeIDs):

        assert kWeibull >= 0.0
        assert mWeibull >= 0.0
        assert minFlawsPerNode > 0
        assert minTotalFlaws > 0 or nodeList.numInternalNodes == 0

        # Prepare the result.
        flaws = VectorDoubleField("Weibull flaw distribution",
                                  nodeList)

        # Construct unique global IDs for all nodes in the NodeList.
        n = numGlobalNodes(nodeList);
        globalIDs = globalNodeIDs(nodeList).internalValues();

        # Only proceed if there are nodes to initialize!
        if n > 0:

            # Compute the Weibull variate parameters.
            beta = (volume/n*kWeibull)**(-1.0/mWeibull)

            # Construct a random number generator.
            g = random.Random(seed)

            # Loop and initialize flaws until:
            # a) every node has the minimum number of flaws per node, and
            # b) we meet the minimum number of total flaws.
            numFlaws = 0;
            numCompletedNodes = 0;
            numFlawsPerNode = [0]*n
            iflaw = 0
            while (numCompletedNodes < n) or (iflaw <= minTotalFlaws):
                iflaw += 1
                if (iflaw % 100 == 0):
                    print("Flaw %i of estimated %i" % (iflaw, n * log(n)))

                # Randomly select a global node.
                iglobal = g.randint(0, n - 1)

                # Increment the number of flaws for this node, and check if this
                # completes this node.
                numFlawsPerNode[iglobal] += 1
                if numFlawsPerNode[iglobal] == minFlawsPerNode:
                    numCompletedNodes += 1;

                # Is this node one of ours?
                if iglobal in globalIDs:

                    i = globalIDs.index(iglobal)
                    assert i >= 0 and i < nodeList.numInternalNodes

                    # The activation energy.
                    epsij = g.weibullvariate(beta, mWeibull)

                    # Add a flaw with this activation energy to this node.
                    flaws[i].append(epsij);

        # Sort the flaws on each node by energy.
        for i in range(nodeList.numInternalNodes):
            v = list(flaws[i])
            v.sort()
            for j in range(len(v)):
                flaws[i][j] = v[j]

        # Post-conditions.
        checkFlaws = 0
        for i in range(nodeList.numInternalNodes):
            nn = len(flaws[i])
            checkFlaws += nn
            assert nn >= minFlawsPerNode
            for j in range(nn - 1):
                assert flaws[i][j] <= flaws[i][j + 1]

        # That's it.
        return flaws

#-------------------------------------------------------------------------------
# 2-D Tensor DamageModel
#-------------------------------------------------------------------------------
class WeibullTensorDamage2d(SolidSpheral.TensorDamageModel2d):

    def __init__(self,
                 nodeList,
                 kWeibull,
                 mWeibull,
                 volume,
                 kernel,
                 seed,
                 strainType,
                 crackGrowthMultiplier = 0.4,
                 minFlawsPerNode = 1,
                 minTotalFlaws = 1):
        self.kWeibull = kWeibull
        self.mWeibull = mWeibull
        self.seed = seed
        SolidSpheral.TensorDamageModel2d.__init__(self,
                                                  nodeList,
                                                  strainType,
                                                  kernel.kernelExtent(),
                                                  crackGrowthMultiplier,
                                                  flawDistribution(seed,
                                                                   kWeibull,
                                                                   mWeibull,
                                                                   volume,
                                                                   nodeList,
                                                                   minFlawsPerNode,
                                                                   minTotalFlaws,
                                                                   SolidSpheral.VectorDoubleField2d,
                                                                   SolidSpheral.numGlobalNodes2d,
                                                                   SolidSpheral.globalNodeIDs2d))
        return

#-------------------------------------------------------------------------------
# 3-D Tensor DamageModel
#-------------------------------------------------------------------------------
class WeibullTensorDamage3d(SolidSpheral.TensorDamageModel3d):

    def __init__(self,
                 nodeList,
                 kWeibull,
                 mWeibull,
                 volume,
                 kernel,
                 seed,
                 strainType,
                 crackGrowthMultiplier = 0.4,
                 minFlawsPerNode = 1,
                 minTotalFlaws = 1):
        self.kWeibull = kWeibull
        self.mWeibull = mWeibull
        self.seed = seed
        SolidSpheral.TensorDamageModel3d.__init__(self,
                                                  nodeList,
                                                  strainType,
                                                  kernel.kernelExtent(),
                                                  crackGrowthMultiplier,
                                                  flawDistribution(seed,
                                                                   kWeibull,
                                                                   mWeibull,
                                                                   volume,
                                                                   nodeList,
                                                                   minFlawsPerNode,
                                                                   minTotalFlaws,
                                                                   SolidSpheral.VectorDoubleField3d,
                                                                   SolidSpheral.numGlobalNodes3d,
                                                                   SolidSpheral.globalNodeIDs3d))
        return

#-------------------------------------------------------------------------------
# 2-D Scalar DamageModel
#-------------------------------------------------------------------------------
class WeibullScalarDamage2d(SolidSpheral.ScalarDamageModel2d):

    def __init__(self,
                 nodeList,
                 damagedNodeList,
                 kWeibull,
                 mWeibull,
                 volume,
                 kernel,
                 seed,
                 crackGrowthMultiplier = 0.4,
                 minFlawsPerNode = 1,
                 minTotalFlaws = 1):
        self.nodeList = nodeList
        self.damagedNodeList = damagedNodeList
        self.kWeibull = kWeibull
        self.mWeibull = mWeibull
        self.seed = seed
        SolidSpheral.ScalarDamageModel2d.__init__(self,
                                                  nodeList,
                                                  damagedNodeList,
                                                  kernel.kernelExtent(),
                                                  crackGrowthMultiplier,
                                                  flawDistribution(seed,
                                                                   kWeibull,
                                                                   mWeibull,
                                                                   volume,
                                                                   nodeList,
                                                                   minFlawsPerNode,
                                                                   minTotalFlaws,
                                                                   SolidSpheral.VectorDoubleField2d,
                                                                   SolidSpheral.numGlobalNodes2d,
                                                                   SolidSpheral.globalNodeIDs2d))
        return

#-------------------------------------------------------------------------------
# 3-D Scalar DamageModel
#-------------------------------------------------------------------------------
class WeibullScalarDamage3d(SolidSpheral.ScalarDamageModel3d):

    def __init__(self,
                 nodeList,
                 damagedNodeList,
                 kWeibull,
                 mWeibull,
                 volume,
                 kernel,
                 seed,
                 crackGrowthMultiplier = 0.4,
                 minFlawsPerNode = 1,
                 minTotalFlaws = 1):
        self.nodeList = nodeList
        self.damagedNodeList = damagedNodeList
        self.kWeibull = kWeibull
        self.mWeibull = mWeibull
        self.seed = seed
        SolidSpheral.ScalarDamageModel3d.__init__(self,
                                                  nodeList,
                                                  damagedNodeList,
                                                  kernel.kernelExtent(),
                                                  crackGrowthMultiplier,
                                                  flawDistribution(seed,
                                                                   kWeibull,
                                                                   mWeibull,
                                                                   volume,
                                                                   nodeList,
                                                                   minFlawsPerNode,
                                                                   minTotalFlaws,
                                                                   SolidSpheral.VectorDoubleField3d,
                                                                   SolidSpheral.numGlobalNodes3d,
                                                                   SolidSpheral.globalNodeIDs3d))
        return

