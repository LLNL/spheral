from Spheral import *
from SpheralTestUtilities import fuzzyEqual

import os
import random

g = random.Random()

WT1d = TableKernel1d(BSplineKernel1d(), 100)
eos1d = GammaLawGasMKS1d(2.0, 2.0)
nodes1d = makeFluidNodeList1d("nodes1d", eos1d)

WT2d = TableKernel2d(BSplineKernel2d(), 100)
eos2d = GammaLawGasMKS2d(2.0, 2.0)
nodes2d = makeFluidNodeList2d("nodes2d", eos2d)

WT3d = TableKernel3d(BSplineKernel3d(), 100)
eos3d = GammaLawGasMKS3d(2.0, 2.0)
nodes3d = makeFluidNodeList3d("nodes3d", eos3d)

n = 10 # 1000
intmin = -2**24
intmax = 2**24
unsignedmin = 0
unsignedmax = 2**32
doublemin = -1e50
doublemax = 1e50
constructor = SidreFileIO

# Size the NodeLists.
nodes1d.numInternalNodes = n
nodes2d.numInternalNodes = n
nodes3d.numInternalNodes = n

#---------------------------------------------------------------------------
    # ScalarFieldListMultipleLists
    #---------------------------------------------------------------------------
fl0 = ScalarFieldList1d()
fl0.copyFields()
fl0.appendNewField("scalar field 1d control", nodes1d, 0.0)
for i in xrange(n):
    fl0[0][i] = g.uniform(doublemin, doublemax)
assert len(fl0) == 1
assert len(fl0[0]) == n

fl1 = ScalarFieldList1d()
fl1.copyFields()
fl1.appendNewField("scalar field 1d control", nodes1d, 0.0)
for i in xrange(n):
    fl1[0][i] = g.uniform(doublemin, doublemax)
assert len(fl1) == 1
assert len(fl1[0]) == n

fl2 = SymTensorFieldList2d()
fl2.copyFields()
fl2.appendNewField("vector field 2d control", nodes2d, SymTensor2d.zero)
for i in xrange(n):
    xx = g.uniform(doublemin, doublemax)
    xy = g.uniform(doublemin, doublemax)
    yy = g.uniform(doublemin, doublemax)
    fl2[0][i] = SymTensor2d(xx, xy,
                            xy, yy)
assert len(fl2) == 1
assert len(fl2[0]) == n

fl3 = VectorFieldList3d()
fl3.copyFields()
fl3.appendNewField("vector field 3d control", nodes3d, Vector3d.zero)
for i in xrange(n):
    fl3[0][i] = Vector3d(g.uniform(doublemin, doublemax),
                         g.uniform(doublemin, doublemax),
                         g.uniform(doublemin, doublemax))
assert len(fl3) == 1
assert len(fl3[0]) == n

fl4 = ScalarFieldList2d()
fl4.copyFields()
fl4.appendNewField("scalar field 2d control", nodes2d, 0.0)
for i in xrange(n):
    fl4[0][i] = g.uniform(doublemin, doublemax)
assert len(fl4) == 1
assert len(fl4[0]) == n

fl5 = IntFieldList1d()
fl5.copyFields()
fl5.appendNewField("int field 1d control", nodes1d, 0)
for i in xrange(n):
    fl5[0][i] = g.randint(intmin, intmax)
assert len(fl5) == 1
assert len(fl5[0]) == n

v0 = ScalarField3d("scalar field 3d control", nodes3d)
for i in xrange(n):
    v0[i] = g.uniform(doublemin, doublemax)
assert len(v0) == n

f = constructor("TestScalarFieldList1d", Write)
f.write(fl0, "FileIOTestBase/ScalarFieldList1d0")
f.write(fl1, "FileIOTestBase/ScalarFieldList1d1")
f.write(fl2, "FileIOTestBase/SymTensorFieldList2d")
f.write(fl3, "FileIOTestBase/VectorFieldList3d")
f.write(fl4, "FileIOTestBase/ScalarFieldList2d")
f.write(fl5, "FileIOTestBase/IntFieldList1d")
f.write(v0, "FileIOTestBase/ScalarField3d")
f.close()

for i in xrange(n):
  print(fl0[0][i]),
print("\n")
for i in xrange(n):
  print(fl1[0][i]),
print("\n")
for i in xrange(n):
  print(fl2[0][i])
print("\n")
for i in xrange(n):
  print(fl3[0][i])
print("\n")
for i in xrange(n):
  print(fl4[0][i]),
print("\n")
for i in xrange(n):
  print(fl5[0][i]),
print("\n")
for i in xrange(n):
  print(v0[i]),
print("\n")
