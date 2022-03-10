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

f = constructor("TestScalarFieldList1d", Read)
fl = ScalarFieldList1d()
fl.copyFields()
fll = ScalarFieldList1d()
fll.copyFields()
fl2r = SymTensorFieldList2d()
fl2r.copyFields()
fl3r = VectorFieldList3d()
fl3r.copyFields()
fl4r = ScalarFieldList2d()
fl4r.copyFields()
fl5r = IntFieldList1d()
fl5r.copyFields()
v = ScalarField3d("scalar field 3d test", nodes3d)

f.read(fl, "FileIOTestBase/ScalarFieldList1d0")
f.read(fll, "FileIOTestBase/ScalarFieldList1d1")
f.read(fl2r, "FileIOTestBase/SymTensorFieldList2d")
f.read(fl3r, "FileIOTestBase/VectorFieldList3d")
f.read(fl4r, "FileIOTestBase/ScalarFieldList2d")
f.read(fl5r, "FileIOTestBase/IntFieldList1d")
f.read(v, "FileIOTestBase/ScalarField3d")
f.close()

# assert len(fl) == len(fl0)
# assert len(fl[0]) == len(fl0[0])
# for i in xrange(n):
#     failUnless(fl[0][i] == fl0[0][i],
#                     "%g != %g @ %i of %i in ScalarFieldList1d0 test" %
#                     (fl[0][i], fl0[0][i], i, n))
# assert len(fll) == len(fl1)
# assert len(fll[0]) == len(fl1[0])
# for i in xrange(n):
#     failUnless(fll[0][i] == fl1[0][i],
#                     "%g != %g @ %i of %i in ScalarFieldList1d1 test" %
#                     (fll[0][i], fl1[0][i], i, n))
# assert len(fl2r) == len(fl2)
# assert len(fl2r[0]) == len(fl2[0])
# for i in xrange(n):
#     failUnless(fl2r[0][i] == fl2[0][i],
#                     "%s != %s @ %i of %i in SymTensorFieldList2d test" %
#                     (str(fl2r[0][i]), str(fl2[0][i]), i, n))
# assert len(fl3r) == len(fl3)
# assert len(fl3r[0]) == len(fl3[0])
# for i in xrange(n):
#     failUnless(fl3r[0][i] == fl3[0][i],
#                     "%s != %s @ %i of %i in VectorFieldList3d test" %
#                     (str(fl3r[0][i]), str(fl3[0][i]), i, n))
# assert len(fl4r) == len(fl4)
# assert len(fl4r[0]) == len(fl4[0])
# for i in xrange(n):
#     failUnless(fl4r[0][i] == fl4[0][i],
#                     "%g != %g @ %i of %i in ScalarFieldList2d test" %
#                     (fl4r[0][i], fl4[0][i], i, n))
# assert len(fl5r) == len(fl5)
# assert len(fl5r[0]) == len(fl5[0])
# for i in xrange(n):
#     failUnless(fl5r[0][i] == fl5[0][i],
#                     "%i != %i @ %i of %i in IntFieldList1d test" %
#                     (fl5r[0][i], fl5[0][i], i, n))
# assert len(v) == len(v0)
# for i in xrange(n):
#     failUnless(v[i] == v0[i],
#                     "%g != %g @ %i of %i in ScalarField3d test" %
#                     (v[i], v0[i], i, n))

for i in xrange(n):
  print(fl[0][i]),
print("\n")
for i in xrange(n):
  print(fll[0][i]),
print("\n")
for i in xrange(n):
  print(fl2r[0][i])
print("\n")
for i in xrange(n):
  print(fl3r[0][i])
print("\n")
for i in xrange(n):
  print(fl4r[0][i]),
print("\n")
for i in xrange(n):
  print(fl5r[0][i]),
print("\n")
for i in xrange(n):
  print(v[i]),
print("\n")
os.remove("TestScalarFieldList1d")
