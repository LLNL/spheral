from Spheral import *

import os
import unittest

def main():
    WT1d = TableKernel1d(BSplineKernel1d(), 100)
    eos1d = GammaLawGasMKS1d(2.0, 2.0)
    nodes1d = makeFluidNodeList1d("nodes1d", eos1d)
    nodes1d.numInternalNodes = 10

    print(ThirdRankTensor1d.nDimensions)
    print(FourthRankTensor1d.nDimensions)
    print(FifthRankTensor1d.nDimensions)

    #x0 = vector_of_int(range(10))
    ##v0 = VectorIntField1d("vector<int> field", nodes1d, x0)
    #v0 = ThirdRankTensorField1d("third rank tensor field 1d control", nodes1d)

    #print(len(v0))

    #print(v0[0].nDimensions)
    #print(v0[1])
    #print(v0[2])
    #print(v0[3])

    #for i in range(3,nodes1d.numInternalNodes):
    #    print(v0[i])




#-------------------------------------------------------------------------------
# Run those tests.
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
    #unittest.main()
