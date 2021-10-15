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

#-------------------------------------------------------------------------------
# Generic FileIO tests.
#-------------------------------------------------------------------------------
class FileIOTestBase:

    #---------------------------------------------------------------------------
    # int
    #---------------------------------------------------------------------------
    def testInt(self):
        x0 = g.randint(self.intmin, self.intmax)
        f = self.constructor("TestInt", Write)
        f.write_int(x0, "FileIOTestBase/a/b/c/d/TestInt")
        f.close()
        f = self.constructor("TestInt", Read)
        x1 = f.read_int("FileIOTestBase/a/b/c/d/TestInt")
        f.close()
        self.failUnless(x1 == x0,
                        "%i != %i in int test" % (x1, x0))
        self.removeFile("TestInt")
        return

    #---------------------------------------------------------------------------
    # bool
    #---------------------------------------------------------------------------
    def testBool(self):
        x0 = g.choice([True, False])
        f = self.constructor("TestBool", Write)
        f.write_bool(x0, "FileIOTestBase/TestBool")
        f.close()
        f = self.constructor("TestBool", Read)
        x1 = f.read_bool("FileIOTestBase/TestBool")
        f.close()
        self.failUnless(x1 == x0,
                        "%s != %s in bool test" % (x1, x0))
        self.removeFile("TestBool")
        return

    #---------------------------------------------------------------------------
    # float
    #---------------------------------------------------------------------------
    def testFloat(self):
        x0 = g.uniform(self.doublemin, self.doublemax)
        f = self.constructor("TestFloat", Write)
        f.write_double(x0, "FileIOTestBase/TestFloat")
        f.close()
        f = self.constructor("TestFloat", Read)
        x1 = f.read_double("FileIOTestBase/TestFloat")
        f.close()
        self.failUnless(x1 == x0,
                        "%s != %s in float test" % (x1, x0))
        self.removeFile("TestFloat")
        return

    #---------------------------------------------------------------------------
    # string
    #---------------------------------------------------------------------------
    def testString(self):
        x0 = "abcdefg"
        f = self.constructor("TestString", Write)
        f.write_string(x0, "FileIOTestBase/TestString")
        f.close()
        f = self.constructor("TestString", Read)
        x1 = f.read_string("FileIOTestBase/TestString")
        f.close()
        self.failUnless(x1 == x0,
                        "%s != %s in string test" % (x1, x0))
        self.removeFile("TestString")
        return

    #---------------------------------------------------------------------------
    # empty string
    #---------------------------------------------------------------------------
    def testEmptyString(self):
        x0 = ""
        f = self.constructor("TestEmptyString", Write)
        f.write_string(x0, "FileIOTestBase/TestEmptyString")
        f.close()
        f = self.constructor("TestEmptyString", Read)
        x1 = f.read_string("FileIOTestBase/TestEmptyString")
        f.close()
        self.failUnless(x1 == x0,
                        "%s != %s in empty string test" % (x1, x0))
        self.removeFile("TestEmptyString")
        return

    #---------------------------------------------------------------------------
    # Vector1d
    #---------------------------------------------------------------------------
    def testVector1d(self):
        x0 = Vector1d(g.uniform(self.doublemin, self.doublemax))
        x1 = Vector1d()
        f = self.constructor("TestVector1d", Write)
        f.write(x0, "FileIOTestBase/TestVector1d")
        f.close()
        f = self.constructor("TestVector1d", Read)
        f.read(x1, "FileIOTestBase/TestVector1d")
        f.close()
        f = self.constructor("TestVector1d-check", Write)
        f.write(x1, "FileIOTestBase/TestVector1d")
        f.close()
        self.failUnless(x1 == x0,
                        "%s != %s in Vector1d test" % (str(x1), str(x0)))
        self.removeFile("TestVector1d")
        self.removeFile("TestVector1d-check")
        return

    #---------------------------------------------------------------------------
    # Tensor1d
    #---------------------------------------------------------------------------
    def testTensor1d(self):
        x0 = Tensor1d(g.uniform(self.doublemin, self.doublemax))
        x1 = Tensor1d()
        f = self.constructor("TestTensor1d", Write)
        f.write(x0, "FileIOTestBase/TestTensor1d")
        f.close()
        f = self.constructor("TestTensor1d", Read)
        f.read(x1, "FileIOTestBase/TestTensor1d")
        f.close()
        self.failUnless(x1 == x0,
                        "%s != %s in Tensor1d test" % (str(x1), str(x0)))
        self.removeFile("TestTensor1d")
        return

    #---------------------------------------------------------------------------
    # SymTensor1d
    #---------------------------------------------------------------------------
    def testSymTensor1d(self):
        x0 = SymTensor1d(g.uniform(self.doublemin, self.doublemax))
        x1 = SymTensor1d()
        f = self.constructor("TestSymTensor1d", Write)
        f.write(x0, "FileIOTestBase/TestSymTensor1d")
        f.close()
        f = self.constructor("TestSymTensor1d", Read)
        f.read(x1, "FileIOTestBase/TestSymTensor1d")
        f.close()
        self.failUnless(x1 == x0,
                        "%s != %s in SymTensor1d test" % (str(x1), str(x0)))
        self.removeFile("TestSymTensor1d")
        return

    #---------------------------------------------------------------------------
    # ThirdRankTensor1d
    #---------------------------------------------------------------------------
    def testThirdRankTensor1d(self):
        x0 = ThirdRankTensor1d()
        for i in xrange(ThirdRankTensor1d.nDimensions):
            for j in xrange(ThirdRankTensor1d.nDimensions):
                for k in xrange(ThirdRankTensor1d.nDimensions):
                    x0(i, j, k, g.uniform(self.doublemin, self.doublemax))
        x1 = ThirdRankTensor1d()
        f = self.constructor("TestThirdRankTensor1d", Write)
        f.write(x0, "FileIOTestBase/TestThirdRankTensor1d")
        f.close()
        f = self.constructor("TestThirdRankTensor1d", Read)
        f.read(x1, "FileIOTestBase/TestThirdRankTensor1d")
        f.close()
        self.failUnless(x1 == x0,
                        "%s != %s in ThirdRankTensor1d test" % (str(x1), str(x0)))
        self.removeFile("TestThirdRankTensor1d")
        return

    #---------------------------------------------------------------------------
    # Vector2d
    #---------------------------------------------------------------------------
    def testVector2d(self):
        x0 = Vector2d(g.uniform(self.doublemin, self.doublemax),
                      g.uniform(self.doublemin, self.doublemax))
        x1 = Vector2d()
        f = self.constructor("TestVector2d", Write)
        f.write(x0, "FileIOTestBase/TestVector2d")
        f.close()
        f = self.constructor("TestVector2d", Read)
        f.read(x1, "FileIOTestBase/TestVector2d")
        f.close()
        self.failUnless(x1 == x0,
                        "%s != %s in Vector2d test" % (str(x1), str(x0)))
        self.removeFile("TestVector2d")
        return

    #---------------------------------------------------------------------------
    # Tensor2d
    #---------------------------------------------------------------------------
    def testTensor2d(self):
        x0 = Tensor2d(g.uniform(self.doublemin, self.doublemax),
                      g.uniform(self.doublemin, self.doublemax),
                      g.uniform(self.doublemin, self.doublemax),
                      g.uniform(self.doublemin, self.doublemax))
        x1 = Tensor2d()
        f = self.constructor("TestTensor2d", Write)
        f.write(x0, "FileIOTestBase/TestTensor2d")
        f.close()
        f = self.constructor("TestTensor2d", Read)
        f.read(x1, "FileIOTestBase/TestTensor2d")
        f.close()
        self.failUnless(x1 == x0,
                        "%s != %s in Tensor2d test" % (str(x1), str(x0)))
        self.removeFile("TestTensor2d")
        return

    #---------------------------------------------------------------------------
    # SymTensor2d
    #---------------------------------------------------------------------------
    def testSymTensor2d(self):
        xy = g.uniform(self.doublemin, self.doublemax)
        x0 = SymTensor2d(g.uniform(self.doublemin, self.doublemax),
                         xy,
                         xy,
                         g.uniform(self.doublemin, self.doublemax))
        x1 = SymTensor2d()
        f = self.constructor("TestSymTensor2d", Write)
        f.write(x0, "FileIOTestBase/TestSymTensor2d")
        f.close()
        f = self.constructor("TestSymTensor2d", Read)
        f.read(x1, "FileIOTestBase/TestSymTensor2d")
        f.close()
        self.failUnless(x1 == x0,
                        "%s != %s in SymTensor2d test" % (str(x1), str(x0)))
        self.removeFile("TestSymTensor2d")
        return

    #---------------------------------------------------------------------------
    # ThirdRankTensor2d
    #---------------------------------------------------------------------------
    def testThirdRankTensor2d(self):
        x0 = ThirdRankTensor2d()
        for i in xrange(ThirdRankTensor2d.nDimensions):
            for j in xrange(ThirdRankTensor2d.nDimensions):
                for k in xrange(ThirdRankTensor2d.nDimensions):
                    x0(i, j, k, g.uniform(self.doublemin, self.doublemax))
        x1 = ThirdRankTensor2d()
        f = self.constructor("TestThirdRankTensor2d", Write)
        f.write(x0, "FileIOTestBase/TestThirdRankTensor2d")
        f.close()
        f = self.constructor("TestThirdRankTensor2d", Read)
        f.read(x1, "FileIOTestBase/TestThirdRankTensor2d")
        f.close()
        self.failUnless(x1 == x0,
                        "%s != %s in ThirdRankTensor2d test" % (str(x1), str(x0)))
        self.removeFile("TestThirdRankTensor2d")
        return

    #---------------------------------------------------------------------------
    # Vector3d
    #---------------------------------------------------------------------------
    def testVector3d(self):
        x0 = Vector3d(g.uniform(self.doublemin, self.doublemax),
                      g.uniform(self.doublemin, self.doublemax),
                      g.uniform(self.doublemin, self.doublemax))
        x1 = Vector3d()
        f = self.constructor("TestVector3d", Write)
        f.write(x0, "FileIOTestBase/TestVector3d")
        f.close()
        f = self.constructor("TestVector3d", Read)
        f.read(x1, "FileIOTestBase/TestVector3d")
        f.close()
        self.failUnless(x1 == x0,
                        "%s != %s in Vector3d test" % (str(x1), str(x0)))
        self.removeFile("TestVector3d")
        return

    #---------------------------------------------------------------------------
    # Tensor3d
    #---------------------------------------------------------------------------
    def testTensor3d(self):
        x0 = Tensor3d(g.uniform(self.doublemin, self.doublemax),
                      g.uniform(self.doublemin, self.doublemax),
                      g.uniform(self.doublemin, self.doublemax),
                      g.uniform(self.doublemin, self.doublemax),
                      g.uniform(self.doublemin, self.doublemax),
                      g.uniform(self.doublemin, self.doublemax),
                      g.uniform(self.doublemin, self.doublemax),
                      g.uniform(self.doublemin, self.doublemax),
                      g.uniform(self.doublemin, self.doublemax))
        x1 = Tensor3d()
        f = self.constructor("TestTensor3d", Write)
        f.write(x0, "FileIOTestBase/TestTensor3d")
        f.close()
        f = self.constructor("TestTensor3d", Read)
        f.read(x1, "FileIOTestBase/TestTensor3d")
        f.close()
        self.failUnless(x1 == x0,
                        "%s != %s in Tensor3d test" % (str(x1), str(x0)))
        self.removeFile("TestTensor3d")
        return

    #---------------------------------------------------------------------------
    # SymTensor3d
    #---------------------------------------------------------------------------
    def testSymTensor3d(self):
        xy = g.uniform(self.doublemin, self.doublemax)
        xz = g.uniform(self.doublemin, self.doublemax)
        yz = g.uniform(self.doublemin, self.doublemax)
        x0 = SymTensor3d(g.uniform(self.doublemin, self.doublemax),
                         xy,
                         xz,
                         xy,
                         g.uniform(self.doublemin, self.doublemax),
                         yz,
                         xz,
                         yz,
                         g.uniform(self.doublemin, self.doublemax))
        x1 = SymTensor3d()
        f = self.constructor("TestSymTensor3d", Write)
        f.write(x0, "FileIOTestBase/TestSymTensor3d")
        f.close()
        f = self.constructor("TestSymTensor3d", Read)
        f.read(x1, "FileIOTestBase/TestSymTensor3d")
        f.close()
        self.failUnless(x1 == x0,
                        "%s != %s in SymTensor3d test" % (str(x1), str(x0)))
        self.removeFile("TestSymTensor3d")
        return

    #---------------------------------------------------------------------------
    # ThirdRankTensor3d
    #---------------------------------------------------------------------------
    def testThirdRankTensor3d(self):
        x0 = ThirdRankTensor3d()
        for i in xrange(ThirdRankTensor3d.nDimensions):
            for j in xrange(ThirdRankTensor3d.nDimensions):
                for k in xrange(ThirdRankTensor3d.nDimensions):
                    x0(i, j, k, g.uniform(self.doublemin, self.doublemax))
        x1 = ThirdRankTensor3d()
        f = self.constructor("TestThirdRankTensor3d", Write)
        f.write(x0, "FileIOTestBase/TestThirdRankTensor3d")
        f.close()
        f = self.constructor("TestThirdRankTensor3d", Read)
        f.read(x1, "FileIOTestBase/TestThirdRankTensor3d")
        f.close()
        self.failUnless(x1 == x0,
                        "%s != %s in ThirdRankTensor3d test" % (str(x1), str(x0)))
        self.removeFile("TestThirdRankTensor3d")
        return

    #---------------------------------------------------------------------------
    # vector<int>
    #---------------------------------------------------------------------------
    def testVectorInt(self):
        for n in (0, self.n):
            filename = "TestVectorInt_%i" % n
            v0 = vector_of_int()
            for i in xrange(n):
                v0.append(g.randint(self.intmin, self.intmax))
            assert len(v0) == n
            f = self.constructor(filename, Write)
            f.write(v0, "FileIOTestBase/vector_of_int")
            f.close()
            f = self.constructor(filename, Read)
            v = vector_of_int()
            f.read(v, "FileIOTestBase/vector_of_int")
            f.close()
            assert len(v) == len(v0)
            for i in xrange(n):
                self.failUnless(v[i] == v0[i],
                                "%i != %i @ %i of %i in vector<int> test" %
                                (v[i], v0[i], i, n))
            self.removeFile(filename)
        return

    #---------------------------------------------------------------------------
    # vector<double>
    #---------------------------------------------------------------------------
    def testVectorDouble(self):
        v0 = vector_of_double()
        for n in (0, self.n):
            filename = "TestVectorDouble_%i" % n
            for i in xrange(n):
                v0.append(g.uniform(self.doublemin, self.doublemax))
            assert len(v0) == n
            f = self.constructor(filename, Write)
            f.write(v0, "FileIOTestBase/vector_of_double")
            f.close()
            f = self.constructor(filename, Read)
            v = vector_of_double()
            f.read(v, "FileIOTestBase/vector_of_double")
            f.close()
            assert len(v) == len(v0)
            for i in xrange(n):
                self.failUnless(v[i] == v0[i],
                                "%g != %g @ %i of %i in vector<double> test" %
                                (v[i], v0[i], i, n))
            self.removeFile(filename)
        return

    #---------------------------------------------------------------------------
    # vector<string>
    #---------------------------------------------------------------------------
    def testVectorString(self):
        chars = "abcdefghijklmnopqrstuvwxyz0123456789!@#$%^&*()-_=+"
        for n in (0, self.n):
            filename = "TestVectorString_%i" % n
            v0 = vector_of_string()
            for i in xrange(n):
                word = ""
                wordlen = g.randint(2, 10)
                for j in xrange(wordlen):
                    word += g.choice(chars)
                v0.append(word)
            assert len(v0) == n
            f = self.constructor(filename, Write)
            f.write(v0, "FileIOTestBase/vector_of_string")
            f.close()
            f = self.constructor(filename, Read)
            v = vector_of_string()
            f.read(v, "FileIOTestBase/vector_of_string")
            f.close()
            assert len(v) == len(v0)
            for i in xrange(n):
                self.failUnless(str(v[i]) == str(v0[i]),
                                ".%s. != .%s. @ %i of %i in vector<string> test" %
                                (v[i], v0[i], i, n))
            self.removeFile(filename)
        return

    #---------------------------------------------------------------------------
    # vector<Vector1d>
    #---------------------------------------------------------------------------
    def testVectorVector1d(self):
        for n in (0, self.n):
            filename = "TestVectorVector1d_%i" % n
            v0 = vector_of_Vector1d()
            for i in xrange(n):
                v0.append(Vector1d(g.uniform(self.doublemin, self.doublemax)))
            assert len(v0) == n
            f = self.constructor(filename, Write)
            f.write(v0, "FileIOTestBase/vector_of_Vector1d")
            f.close()
            f = self.constructor(filename, Read)
            v = vector_of_Vector1d()
            f.read(v, "FileIOTestBase/vector_of_Vector1d")
            f.close()
            assert len(v) == len(v0)
            for i in xrange(n):
                self.failUnless(v[i] == v0[i],
                                "%s != %s @ %i of %i in vector<Vector1d> test" %
                                (str(v[i]), str(v0[i]), i, n))
            self.removeFile(filename)
        return

    #---------------------------------------------------------------------------
    # vector<Vector2d>
    #---------------------------------------------------------------------------
    def testVectorVector2d(self):
        for n in (0, self.n):
            filename = "TestVectorVector2d_%i" % n
            v0 = vector_of_Vector2d()
            for i in xrange(n):
                v0.append(Vector2d(g.uniform(self.doublemin, self.doublemax),
                                   g.uniform(self.doublemin, self.doublemax)))
            assert len(v0) == n
            f = self.constructor(filename, Write)
            f.write(v0, "FileIOTestBase/vector_of_Vector2d")
            f.close()
            f = self.constructor(filename, Read)
            v1 = vector_of_Vector2d()
            f.read(v1, "FileIOTestBase/vector_of_Vector2d")
            f.close()
            assert len(v1) == len(v0)
            for i in xrange(n):
                self.failUnless(fuzzyEqual(list(v1[i]), list(v0[i])),
                                "%s != %s @ %i of %i in vector<Vector2d> test" %
                                (str(v1[i]), str(v0[i]), i, n))
            self.removeFile(filename)
        return

    #---------------------------------------------------------------------------
    # Vector<Vector3d>
    #---------------------------------------------------------------------------
    def testVectorVector3d(self):
        for n in (0, self.n):
            filename = "TestVectorVector3d_%i" % n
            v0 = vector_of_Vector3d()
            for i in xrange(n):
                v0.append(Vector3d(g.uniform(self.doublemin, self.doublemax),
                                   g.uniform(self.doublemin, self.doublemax),
                                   g.uniform(self.doublemin, self.doublemax)))
            assert len(v0) == n
            f = self.constructor(filename, Write)
            f.write(v0, "FileIOTestBase/vector_of_Vector3d")
            f.close()
            f = self.constructor(filename, Read)
            v = vector_of_Vector3d()
            f.read(v, "FileIOTestBase/vector_of_Vector3d")
            f.close()
            assert len(v) == len(v0)
            for i in xrange(n):
                self.failUnless(fuzzyEqual(list(v[i]), list(v0[i])),
                                "%s != %s @ %i of %i in vector<Vector3d> test" %
                                (str(v[i]), str(v0[i]), i, n))
            self.removeFile(filename)
        return

    #---------------------------------------------------------------------------
    # vector<Tensor1d>
    #---------------------------------------------------------------------------
    def testVectorTensor1d(self):
        for n in (0, self.n):
            filename = "TestVectorTensor1d_%i" % n
            v0 = vector_of_Tensor1d()
            for i in xrange(n):
                v0.append(Tensor1d(g.uniform(self.doublemin, self.doublemax)))
            assert len(v0) == n
            f = self.constructor(filename, Write)
            f.write(v0, "FileIOTestBase/vector_of_Tensor1d")
            f.close()
            f = self.constructor(filename, Read)
            v = vector_of_Tensor1d()
            f.read(v, "FileIOTestBase/vector_of_Tensor1d")
            f.close()
            assert len(v) == len(v0)
            for i in xrange(n):
                self.failUnless(fuzzyEqual(list(v[i]), list(v0[i])),
                                "%s != %s @ %i of %i in vector<Tensor1d> test" %
                                (str(v[i]), str(v0[i]), i, n))
            self.removeFile(filename)
        return

    #---------------------------------------------------------------------------
    # vector<Tensor2d>
    #---------------------------------------------------------------------------
    def testVectorTensor2d(self):
        for n in (0, self.n):
            filename = "TestVectorTensor2d_%i" % n
            v0 = vector_of_Tensor2d()
            for i in xrange(n):
                v0.append(Tensor2d(g.uniform(self.doublemin, self.doublemax),
                                   g.uniform(self.doublemin, self.doublemax),
                                   g.uniform(self.doublemin, self.doublemax),
                                   g.uniform(self.doublemin, self.doublemax)))
            assert len(v0) == n
            f = self.constructor(filename, Write)
            f.write(v0, "FileIOTestBase/vector_of_Tensor2d")
            f.close()
            f = self.constructor(filename, Read)
            v = vector_of_Tensor2d()
            f.read(v, "FileIOTestBase/vector_of_Tensor2d")
            f.close()
            assert len(v) == len(v0)
            for i in xrange(n):
                self.failUnless(fuzzyEqual(list(v[i]), list(v0[i])),
                                "%s != %s @ %i of %i in vector<Tensor2d> test" %
                                (str(v[i]), str(v0[i]), i, n))
            self.removeFile(filename)
        return

    #---------------------------------------------------------------------------
    # vector<Tensor3d>
    #---------------------------------------------------------------------------
    def testVectorTensor3d(self):
        for n in (0, self.n):
            filename = "TestVectorTensor3d_%i" % n
            v0 = vector_of_Tensor3d()
            for i in xrange(n):
                v0.append(Tensor3d(g.uniform(self.doublemin, self.doublemax),
                                   g.uniform(self.doublemin, self.doublemax),
                                   g.uniform(self.doublemin, self.doublemax),
                                   g.uniform(self.doublemin, self.doublemax),
                                   g.uniform(self.doublemin, self.doublemax),
                                   g.uniform(self.doublemin, self.doublemax),
                                   g.uniform(self.doublemin, self.doublemax),
                                   g.uniform(self.doublemin, self.doublemax),
                                   g.uniform(self.doublemin, self.doublemax)))
            assert len(v0) == n
            f = self.constructor(filename, Write)
            f.write(v0, "FileIOTestBase/vector_of_Tensor3d")
            f.close()
            f = self.constructor(filename, Read)
            v = vector_of_Tensor3d()
            f.read(v, "FileIOTestBase/vector_of_Tensor3d")
            f.close()
            assert len(v) == len(v0)
            for i in xrange(n):
                self.failUnless(fuzzyEqual(list(v[i]), list(v0[i])),
                                "%s != %s @ %i of %i in vector<Tensor3d> test" %
                                (str(v[i]), str(v0[i]), i, n))
            self.removeFile(filename)
        return

    #---------------------------------------------------------------------------
    # vector<SymTensor1d>
    #---------------------------------------------------------------------------
    def testVectorSymTensor1d(self):
        for n in (0, self.n):
            filename = "TestVectorSymTensor1d_%i" % n
            v0 = vector_of_SymTensor1d()
            for i in xrange(n):
                v0.append(SymTensor1d(g.uniform(self.doublemin, self.doublemax)))
            assert len(v0) == n
            f = self.constructor(filename, Write)
            f.write(v0, "FileIOTestBase/vector_of_SymTensor1d")
            f.close()
            f = self.constructor(filename, Read)
            v = vector_of_SymTensor1d()
            f.read(v, "FileIOTestBase/vector_of_SymTensor1d")
            f.close()
            assert len(v) == len(v0)
            for i in xrange(n):
                self.failUnless(fuzzyEqual(list(v[i]), list(v0[i])),
                                "%s != %s @ %i of %i in vector<SymTensor1d> test" %
                                (str(v[i]), str(v0[i]), i, n))
            self.removeFile(filename)
        return

    #---------------------------------------------------------------------------
    # vector<SymTensor2d>
    #---------------------------------------------------------------------------
    def testVectorSymTensor2d(self):
        for n in (0, self.n):
            filename = "TestVectorSymTensor2d_%i" % n
            v0 = vector_of_SymTensor2d()
            for i in xrange(n):
                xx = g.uniform(self.doublemin, self.doublemax)
                xy = g.uniform(self.doublemin, self.doublemax)
                yy = g.uniform(self.doublemin, self.doublemax)
                v0.append(SymTensor2d(xx, xy,
                                      xy, yy))
            assert len(v0) == n
            f = self.constructor(filename, Write)
            f.write(v0, "FileIOTestBase/vector_of_SymTensor2d")
            f.close()
            f = self.constructor(filename, Read)
            v = vector_of_SymTensor2d()
            f.read(v, "FileIOTestBase/vector_of_SymTensor2d")
            f.close()
            assert len(v) == len(v0)
            for i in xrange(n):
                self.failUnless(fuzzyEqual(list(v[i]), list(v0[i])),
                                "%s != %s @ %i of %i in vector<SymTensor2d> test" %
                                (str(v[i]), str(v0[i]), i, n))
            self.removeFile(filename)
        return

    #---------------------------------------------------------------------------
    # vector<SymTensor3d>
    #---------------------------------------------------------------------------
    def testVectorSymTensor3d(self):
        for n in (0, self.n):
            filename = "TestVectorSymTensor3d_%i" % n
            v0 = vector_of_SymTensor3d()
            for i in xrange(n):
                xx = g.uniform(self.doublemin, self.doublemax)
                xy = g.uniform(self.doublemin, self.doublemax)
                xz = g.uniform(self.doublemin, self.doublemax)
                yy = g.uniform(self.doublemin, self.doublemax)
                yz = g.uniform(self.doublemin, self.doublemax)
                zz = g.uniform(self.doublemin, self.doublemax)
                v0.append(SymTensor3d(xx, xy, xz,
                                      xy, yy, yz,
                                      xz, yz, zz))
            assert len(v0) == n
            f = self.constructor(filename, Write)
            f.write(v0, "FileIOTestBase/vector_of_SymTensor3d")
            f.close()
            f = self.constructor(filename, Read)
            v = vector_of_SymTensor3d()
            f.read(v, "FileIOTestBase/vector_of_SymTensor3d")
            f.close()
            assert len(v) == len(v0)
            for i in xrange(n):
                self.failUnless(fuzzyEqual(list(v[i]), list(v0[i])),
                                "%s != %s @ %i of %i in vector<SymTensor3d> test" %
                                (str(v[i]), str(v0[i]), i, n))
            self.removeFile(filename)
        return

    #---------------------------------------------------------------------------
    # vector<ThirdRankTensor1d>
    #---------------------------------------------------------------------------
    def testVectorThirdRankTensor1d(self):
        for n in (0, self.n):
            filename = "TestVectorThirdRankTensor1d_%i" % n
            v0 = vector_of_ThirdRankTensor1d()
            for i in xrange(n):
                val = ThirdRankTensor1d()
                for i in xrange(ThirdRankTensor1d.nDimensions):
                    for j in xrange(ThirdRankTensor1d.nDimensions):
                        for k in xrange(ThirdRankTensor1d.nDimensions):
                            val(i, j, k, g.uniform(self.doublemin, self.doublemax))
                v0.append(val)
            assert len(v0) == n
            f = self.constructor(filename, Write)
            f.write(v0, "FileIOTestBase/vector_of_ThirdRankTensor1d")
            f.close()
            f = self.constructor(filename, Read)
            v = vector_of_ThirdRankTensor1d()
            f.read(v, "FileIOTestBase/vector_of_ThirdRankTensor1d")
            f.close()
            assert len(v) == len(v0)
            for i in xrange(n):
                self.failUnless(fuzzyEqual(list(v[i]), list(v0[i])),
                                "%s != %s @ %i of %i in vector<ThirdRankTensor1d> test" %
                                (str(v[i]), str(v0[i]), i, n))
            self.removeFile(filename)
        return

    #---------------------------------------------------------------------------
    # vector<ThirdRankTensor2d>
    #---------------------------------------------------------------------------
    def testVectorThirdRankTensor2d(self):
        for n in (0, self.n):
            filename = "TestVectorThirdRankTensor2d_%i" % n
            v0 = vector_of_ThirdRankTensor2d()
            for i in xrange(n):
                val = ThirdRankTensor2d()
                for i in xrange(ThirdRankTensor2d.nDimensions):
                    for j in xrange(ThirdRankTensor2d.nDimensions):
                        for k in xrange(ThirdRankTensor2d.nDimensions):
                            val(i, j, k, g.uniform(self.doublemin, self.doublemax))
                v0.append(val)
            assert len(v0) == n
            f = self.constructor(filename, Write)
            f.write(v0, "FileIOTestBase/vector_of_ThirdRankTensor2d")
            f.close()
            f = self.constructor(filename, Read)
            v = vector_of_ThirdRankTensor2d()
            f.read(v, "FileIOTestBase/vector_of_ThirdRankTensor2d")
            f.close()
            assert len(v) == len(v0)
            for i in xrange(n):
                self.failUnless(fuzzyEqual(list(v[i]), list(v0[i])),
                                "%s != %s @ %i of %i in vector<ThirdRankTensor2d> test" %
                                (str(v[i]), str(v0[i]), i, n))
            self.removeFile(filename)
        return

    #---------------------------------------------------------------------------
    # vector<ThirdRankTensor3d>
    #---------------------------------------------------------------------------
    def testVectorThirdRankTensor3d(self):
        for n in (0, self.n):
            filename = "TestVectorThirdRankTensor3d_%i" % n
            v0 = vector_of_ThirdRankTensor3d()
            for i in xrange(n):
                val = ThirdRankTensor3d()
                for i in xrange(ThirdRankTensor2d.nDimensions):
                    for j in xrange(ThirdRankTensor2d.nDimensions):
                        for k in xrange(ThirdRankTensor2d.nDimensions):
                            val(i, j, k, g.uniform(self.doublemin, self.doublemax))
                v0.append(val)
            assert len(v0) == n
            f = self.constructor(filename, Write)
            f.write(v0, "FileIOTestBase/vector_of_ThirdRankTensor3d")
            f.close()
            f = self.constructor(filename, Read)
            v = vector_of_ThirdRankTensor3d()
            f.read(v, "FileIOTestBase/vector_of_ThirdRankTensor3d")
            f.close()
            assert len(v) == len(v0)
            for i in xrange(n):
                self.failUnless(fuzzyEqual(list(v[i]), list(v0[i])),
                                "%s != %s @ %i of %i in vector<ThirdRankTensor3d> test" %
                                (str(v[i]), str(v0[i]), i, n))
            self.removeFile(filename)
        return

    #---------------------------------------------------------------------------
    # Intfield1d
    #---------------------------------------------------------------------------
    def testIntField1d(self):
        v0 = IntField1d("int field 1d control", nodes1d)
        for i in xrange(self.n):
            v0[i] = g.randint(self.intmin, self.intmax)
        assert len(v0) == self.n
        f = self.constructor("TestIntField1d", Write)
        f.write(v0, "FileIOTestBase/IntField1d")
        f.close()
        f = self.constructor("TestIntField1d", Read)
        v = IntField1d("int field 1d test", nodes1d)
        f.read(v, "FileIOTestBase/IntField1d")
        f.close()
        assert len(v) == len(v0)
        for i in xrange(self.n):
            self.failUnless(v[i] == v0[i],
                            "%i != %i @ %i of %i in IntField1d test" %
                            (v[i], v0[i], i, self.n))
        self.removeFile("TestIntField1d")
        return

    #---------------------------------------------------------------------------
    # IntField2d
    #---------------------------------------------------------------------------
    def testIntField2d(self):
        v0 = IntField2d("int field 2d control", nodes2d)
        for i in xrange(self.n):
            v0[i] = g.randint(self.intmin, self.intmax)
        assert len(v0) == self.n
        f = self.constructor("TestIntField2d", Write)
        f.write(v0, "FileIOTestBase/IntField2d")
        f.close()
        f = self.constructor("TestIntField2d", Read)
        v = IntField2d("int field 2d test", nodes2d)
        f.read(v, "FileIOTestBase/IntField2d")
        f.close()
        assert len(v) == len(v0)
        for i in xrange(self.n):
            self.failUnless(v[i] == v0[i],
                            "%i != %i @ %i of %i in IntField2d test" %
                            (v[i], v0[i], i, self.n))
        self.removeFile("TestIntField2d")
        return

    #---------------------------------------------------------------------------
    # IntField3d
    #---------------------------------------------------------------------------
    def testIntField3d(self):
        v0 = IntField3d("int field 3d control", nodes3d)
        for i in xrange(self.n):
            v0[i] = g.randint(self.intmin, self.intmax)
        assert len(v0) == self.n
        f = self.constructor("TestIntField3d", Write)
        f.write(v0, "FileIOTestBase/IntField3d")
        f.close()
        f = self.constructor("TestIntField3d", Read)
        v = IntField3d("int field 3d test", nodes3d)
        f.read(v, "FileIOTestBase/IntField3d")
        f.close()
        assert len(v) == len(v0)
        for i in xrange(self.n):
            self.failUnless(v[i] == v0[i],
                            "%i != %i @ %i of %i in IntField3d test" %
                            (v[i], v0[i], i, self.n))
        self.removeFile("TestIntField3d")
        return

    #---------------------------------------------------------------------------
    # UnsignedField1d
    #---------------------------------------------------------------------------
    def testUnsignedField1d(self):
        v0 = UnsignedField1d("unsigned field 1d control", nodes1d)
        for i in xrange(self.n):
            v0[i] = g.randint(self.unsignedmin, self.unsignedmax)
        assert len(v0) == self.n
        f = self.constructor("TestUnsignedField1d", Write)
        f.write(v0, "FileIOTestBase/UnsignedField1d")
        f.close()
        f = self.constructor("TestUnsignedField1d", Read)
        v = UnsignedField1d("unsigned field 1d test", nodes1d)
        f.read(v, "FileIOTestBase/UnsignedField1d")
        f.close()
        assert len(v) == len(v0)
        for i in xrange(self.n):
            self.failUnless(v[i] == v0[i],
                            "%i != %i @ %i of %i in UnsignedField1d test" %
                            (v[i], v0[i], i, self.n))
        self.removeFile("TestUnsignedField1d")
        return

    #---------------------------------------------------------------------------
    # ScalarField1d
    #---------------------------------------------------------------------------
    def testScalarField1d(self):
        v0 = ScalarField1d("scalar field 1d control", nodes1d)
        for i in xrange(self.n):
            v0[i] = g.uniform(self.doublemin, self.doublemax)
        assert len(v0) == self.n
        f = self.constructor("TestScalarField1d", Write)
        f.write(v0, "FileIOTestBase/ScalarField1d")
        f.close()
        f = self.constructor("TestScalarField1d", Read)
        v = ScalarField1d("scalar field 1d test", nodes1d)
        f.read(v, "FileIOTestBase/ScalarField1d")
        f.close()
        assert len(v) == len(v0)
        for i in xrange(self.n):
            self.failUnless(v[i] == v0[i],
                            "%g != %g @ %i of %i in ScalarField1d test" %
                            (v[i], v0[i], i, self.n))
        self.removeFile("TestScalarField1d")
        return

    #---------------------------------------------------------------------------
    # ScalarField2d
    #---------------------------------------------------------------------------
    def testScalarField2d(self):
        v0 = ScalarField2d("scalar field 2d control", nodes2d)
        for i in xrange(self.n):
            v0[i] = g.uniform(self.doublemin, self.doublemax)
        assert len(v0) == self.n
        f = self.constructor("TestScalarField2d", Write)
        f.write(v0, "FileIOTestBase/ScalarField2d")
        f.close()
        f = self.constructor("TestScalarField2d", Read)
        v = ScalarField2d("scalar field 2d test", nodes2d)
        f.read(v, "FileIOTestBase/ScalarField2d")
        f.close()
        assert len(v) == len(v0)
        for i in xrange(self.n):
            self.failUnless(v[i] == v0[i],
                            "%g != %g @ %i of %i in ScalarField2d test" %
                            (v[i], v0[i], i, self.n))
        self.removeFile("TestScalarField2d")
        return

    #---------------------------------------------------------------------------
    # ScalarField3d
    #---------------------------------------------------------------------------
    def testScalarField3d(self):
        v0 = ScalarField3d("scalar field 3d control", nodes3d)
        for i in xrange(self.n):
            v0[i] = g.uniform(self.doublemin, self.doublemax)
        assert len(v0) == self.n
        f = self.constructor("TestScalarField3d", Write)
        f.write(v0, "FileIOTestBase/ScalarField3d")
        f.close()
        f = self.constructor("TestScalarField3d", Read)
        v = ScalarField3d("scalar field 3d test", nodes3d)
        f.read(v, "FileIOTestBase/ScalarField3d")
        f.close()
        assert len(v) == len(v0)
        for i in xrange(self.n):
            self.failUnless(v[i] == v0[i],
                            "%g != %g @ %i of %i in ScalarField3d test" %
                            (v[i], v0[i], i, self.n))
        self.removeFile("TestScalarField3d")
        return

    #---------------------------------------------------------------------------
    # VectorField1d
    #---------------------------------------------------------------------------
    def testVectorField1d(self):
        v0 = VectorField1d("vector field 1d control", nodes1d)
        for i in xrange(self.n):
            v0[i] = Vector1d(g.uniform(self.doublemin, self.doublemax))
        assert len(v0) == self.n
        f = self.constructor("TestVectorField1d", Write)
        f.write(v0, "FileIOTestBase/VectorField1d")
        f.close()
        f = self.constructor("TestVectorField1d", Read)
        v = VectorField1d("vector field 1d test", nodes1d)
        f.read(v, "FileIOTestBase/VectorField1d")
        f.close()
        assert len(v) == len(v0)
        for i in xrange(self.n):
            self.failUnless(v[i] == v0[i],
                            "%s != %s @ %i of %i in VectorField1d test" %
                            (str(v[i]), str(v0[i]), i, self.n))
        self.removeFile("TestVectorField1d")
        return

    #---------------------------------------------------------------------------
    # VectorField2d
    #---------------------------------------------------------------------------
    def testVectorField2d(self):
        v0 = VectorField2d("vector field 2d control", nodes2d)
        for i in xrange(self.n):
            v0[i] = Vector2d(g.uniform(self.doublemin, self.doublemax),
                             g.uniform(self.doublemin, self.doublemax))
        assert len(v0) == self.n
        f = self.constructor("TestVectorField2d", Write)
        f.write(v0, "FileIOTestBase/VectorField2d")
        f.close()
        f = self.constructor("TestVectorField2d", Read)
        v = VectorField2d("vector field 2d test", nodes2d)
        f.read(v, "FileIOTestBase/VectorField2d")
        f.close()
        assert len(v) == len(v0)
        for i in xrange(self.n):
            self.failUnless(v[i] == v0[i],
                            "%s != %s @ %i of %i in VectorField2d test" %
                            (str(v[i]), str(v0[i]), i, self.n))
        self.removeFile("TestVectorField2d")
        return

    #---------------------------------------------------------------------------
    # VectorField3d
    #---------------------------------------------------------------------------
    def testVectorField3d(self):
        v0 = VectorField3d("vector field 3d control", nodes3d)
        for i in xrange(self.n):
            v0[i] = Vector3d(g.uniform(self.doublemin, self.doublemax),
                             g.uniform(self.doublemin, self.doublemax),
                             g.uniform(self.doublemin, self.doublemax))
        assert len(v0) == self.n
        f = self.constructor("TestVectorField3d", Write)
        f.write(v0, "FileIOTestBase/VectorField3d")
        f.close()
        f = self.constructor("TestVectorField3d", Read)
        v = VectorField3d("vector field 3d test", nodes3d)
        f.read(v, "FileIOTestBase/VectorField3d")
        f.close()
        assert len(v) == len(v0)
        for i in xrange(self.n):
            self.failUnless(v[i] == v0[i],
                            "%s != %s @ %i of %i in VectorField3d test" %
                            (str(v[i]), str(v0[i]), i, self.n))
        self.removeFile("TestVectorField3d")
        return

    #---------------------------------------------------------------------------
    # TensorField1d
    #---------------------------------------------------------------------------
    def testTensorField1d(self):
        v0 = TensorField1d("tensor field 1d control", nodes1d)
        for i in xrange(self.n):
            v0[i] = Tensor1d(g.uniform(self.doublemin, self.doublemax))
        assert len(v0) == self.n
        f = self.constructor("TestTensorField1d", Write)
        f.write(v0, "FileIOTestBase/TensorField1d")
        f.close()
        f = self.constructor("TestTensorField1d", Read)
        v = TensorField1d("tensor field 1d test", nodes1d)
        f.read(v, "FileIOTestBase/TensorField1d")
        f.close()
        assert len(v) == len(v0)
        for i in xrange(self.n):
            self.failUnless(v[i] == v0[i],
                            "%s != %s @ %i of %i in TensorField1d test" %
                            (str(v[i]), str(v0[i]), i, self.n))
        self.removeFile("TestTensorField1d")
        return

    #---------------------------------------------------------------------------
    # TensorField2d
    #---------------------------------------------------------------------------
    def testTensorField2d(self):
        v0 = TensorField2d("tensor field 2d control", nodes2d)
        for i in xrange(self.n):
            v0[i] = Tensor2d(g.uniform(self.doublemin, self.doublemax),
                             g.uniform(self.doublemin, self.doublemax),
                             g.uniform(self.doublemin, self.doublemax),
                             g.uniform(self.doublemin, self.doublemax))
        assert len(v0) == self.n
        f = self.constructor("TestTensorField2d", Write)
        f.write(v0, "FileIOTestBase/TensorField2d")
        f.close()
        f = self.constructor("TestTensorField2d", Read)
        v = TensorField2d("tensor field 2d test", nodes2d)
        f.read(v, "FileIOTestBase/TensorField2d")
        f.close()
        assert len(v) == len(v0)
        for i in xrange(self.n):
            self.failUnless(v[i] == v0[i],
                            "%s != %s @ %i of %i in TensorField2d test" %
                            (str(v[i]), str(v0[i]), i, self.n))
        self.removeFile("TestTensorField2d")
        return

    #---------------------------------------------------------------------------
    # TensorField3d
    #---------------------------------------------------------------------------
    def testTensorField3d(self):
        v0 = TensorField3d("tensor field 3d control", nodes3d)
        for i in xrange(self.n):
            v0[i] = Tensor3d(g.uniform(self.doublemin, self.doublemax),
                             g.uniform(self.doublemin, self.doublemax),
                             g.uniform(self.doublemin, self.doublemax),
                             g.uniform(self.doublemin, self.doublemax),
                             g.uniform(self.doublemin, self.doublemax),
                             g.uniform(self.doublemin, self.doublemax),
                             g.uniform(self.doublemin, self.doublemax),
                             g.uniform(self.doublemin, self.doublemax),
                             g.uniform(self.doublemin, self.doublemax))
        assert len(v0) == self.n
        f = self.constructor("TestTensorField3d", Write)
        f.write(v0, "FileIOTestBase/TensorField3d")
        f.close()
        f = self.constructor("TestTensorField3d", Read)
        v = TensorField3d("tensor field 3d test", nodes3d)
        f.read(v, "FileIOTestBase/TensorField3d")
        f.close()
        assert len(v) == len(v0)
        for i in xrange(self.n):
            self.failUnless(v[i] == v0[i],
                            "%s != %s @ %i of %i in TensorField3d test" %
                            (str(v[i]), str(v0[i]), i, self.n))
        self.removeFile("TestTensorField3d")
        return

    #---------------------------------------------------------------------------
    # SymTensorField1d
    #---------------------------------------------------------------------------
    def testSymTensorField1d(self):
        v0 = SymTensorField1d("symtensor field 1d control", nodes1d)
        for i in xrange(self.n):
            v0[i] = SymTensor1d(g.uniform(self.doublemin, self.doublemax))
        assert len(v0) == self.n
        f = self.constructor("TestSymTensorField1d", Write)
        f.write(v0, "FileIOTestBase/SymTensorField1d")
        f.close()
        f = self.constructor("TestSymTensorField1d", Read)
        v = SymTensorField1d("symtensor field 1d test", nodes1d)
        f.read(v, "FileIOTestBase/SymTensorField1d")
        f.close()
        assert len(v) == len(v0)
        for i in xrange(self.n):
            self.failUnless(v[i] == v0[i],
                            "%s != %s @ %i of %i in SymTensorField1d test" %
                            (str(v[i]), str(v0[i]), i, self.n))
        self.removeFile("TestSymTensorField1d")
        return

    #---------------------------------------------------------------------------
    # SymTensorField2d
    #---------------------------------------------------------------------------
    def testSymTensorField2d(self):
        v0 = SymTensorField2d("symtensor field 2d control", nodes2d)
        for i in xrange(self.n):
            xx = g.uniform(self.doublemin, self.doublemax)
            xy = g.uniform(self.doublemin, self.doublemax)
            yy = g.uniform(self.doublemin, self.doublemax)
            v0[i] = SymTensor2d(xx, xy,
                                xy, yy)
        assert len(v0) == self.n
        f = self.constructor("TestSymTensorField2d", Write)
        f.write(v0, "FileIOTestBase/SymTensorField2d")
        f.close()
        f = self.constructor("TestSymTensorField2d", Read)
        v = SymTensorField2d("symtensor field 2d test", nodes2d)
        f.read(v, "FileIOTestBase/SymTensorField2d")
        f.close()
        assert len(v) == len(v0)
        for i in xrange(self.n):
            self.failUnless(v[i] == v0[i],
                            "%s != %s @ %i of %i in SymTensorField2d test" %
                            (str(v[i]), str(v0[i]), i, self.n))
        self.removeFile("TestSymTensorField2d")
        return

    #---------------------------------------------------------------------------
    # SymTensorField3d
    #---------------------------------------------------------------------------
    def testSymTensorField3d(self):
        v0 = SymTensorField3d("symtensor field 3d control", nodes3d)
        for i in xrange(self.n):
            xx = g.uniform(self.doublemin, self.doublemax)
            xy = g.uniform(self.doublemin, self.doublemax)
            xz = g.uniform(self.doublemin, self.doublemax)
            yy = g.uniform(self.doublemin, self.doublemax)
            yz = g.uniform(self.doublemin, self.doublemax)
            zz = g.uniform(self.doublemin, self.doublemax)
            v0[i] = SymTensor3d(xx, xy, xz,
                                xy, yy, yz,
                                xz, yz, zz)
        assert len(v0) == self.n
        f = self.constructor("TestSymTensorField3d", Write)
        f.write(v0, "FileIOTestBase/SymTensorField3d")
        f.close()
        f = self.constructor("TestSymTensorField3d", Read)
        v = SymTensorField3d("symtensor field 3d test", nodes3d)
        f.read(v, "FileIOTestBase/SymTensorField3d")
        f.close()
        assert len(v) == len(v0)
        for i in xrange(self.n):
            self.failUnless(v[i] == v0[i],
                            "%s != %s @ %i of %i in SymTensorField3d test" %
                            (str(v[i]), str(v0[i]), i, self.n))
        self.removeFile("TestSymTensorField3d")
        return

    #---------------------------------------------------------------------------
    # ThirdRankTensorField1d
    #---------------------------------------------------------------------------
    def testThirdRankTensorField1d(self):
        v0 = ThirdRankTensorField1d("third rank tensor field 1d control", nodes1d)
        for ii in xrange(self.n):
            v0[ii] = ThirdRankTensor1d()
            for i in xrange(ThirdRankTensor1d.nDimensions):
                for j in xrange(ThirdRankTensor1d.nDimensions):
                    for k in xrange(ThirdRankTensor1d.nDimensions):
                        v0[ii](i, j, k, g.uniform(self.doublemin, self.doublemax))
        assert len(v0) == self.n
        f = self.constructor("TestThirdRankTensorField1d", Write)
        f.write(v0, "FileIOTestBase/ThirdRankTensorField1d")
        f.close()
        f = self.constructor("TestThirdRankTensorField1d", Read)
        v = ThirdRankTensorField1d("third rank tensor field 1d test", nodes1d)
        f.read(v, "FileIOTestBase/ThirdRankTensorField1d")
        f.close()
        assert len(v) == len(v0)
        for i in xrange(self.n):
            self.failUnless(v[i] == v0[i],
                            "%s != %s @ %i of %i in ThirdRankTensorField1d test" %
                            (str(v[i]), str(v0[i]), i, self.n))
        self.removeFile("TestThirdRankTensorField1d")
        return

    #---------------------------------------------------------------------------
    # ThirdRankTensorField2d
    #---------------------------------------------------------------------------
    def testThirdRankTensorField2d(self):
        v0 = ThirdRankTensorField2d("third rank tensor field 2d control", nodes2d)
        for ii in xrange(self.n):
            v0[ii] = ThirdRankTensor2d()
            for i in xrange(ThirdRankTensor2d.nDimensions):
                for j in xrange(ThirdRankTensor2d.nDimensions):
                    for k in xrange(ThirdRankTensor2d.nDimensions):
                        v0[ii](i, j, k, g.uniform(self.doublemin, self.doublemax))
        assert len(v0) == self.n
        f = self.constructor("TestThirdRankTensorField2d", Write)
        f.write(v0, "FileIOTestBase/ThirdRankTensorField2d")
        f.close()
        f = self.constructor("TestThirdRankTensorField2d", Read)
        v = ThirdRankTensorField2d("third rank tensor field 2d test", nodes2d)
        f.read(v, "FileIOTestBase/ThirdRankTensorField2d")
        f.close()
        assert len(v) == len(v0)
        for i in xrange(self.n):
            self.failUnless(v[i] == v0[i],
                            "%s != %s @ %i of %i in ThirdRankTensorField2d test" %
                            (str(v[i]), str(v0[i]), i, self.n))
        self.removeFile("TestThirdRankTensorField2d")
        return

    #---------------------------------------------------------------------------
    # ThirdRankTensorField3d
    #---------------------------------------------------------------------------
    def testThirdRankTensorField3d(self):
        v0 = ThirdRankTensorField3d("third rank tensor field 3d control", nodes3d)
        for ii in xrange(self.n):
            v0[ii] = ThirdRankTensor3d()
            for i in xrange(ThirdRankTensor3d.nDimensions):
                for j in xrange(ThirdRankTensor3d.nDimensions):
                    for k in xrange(ThirdRankTensor3d.nDimensions):
                        v0[ii](i, j, k, g.uniform(self.doublemin, self.doublemax))
        assert len(v0) == self.n
        f = self.constructor("TestThirdRankTensorField3d", Write)
        f.write(v0, "FileIOTestBase/ThirdRankTensorField3d")
        f.close()
        f = self.constructor("TestThirdRankTensorField3d", Read)
        v = ThirdRankTensorField3d("third rank tensor field 3d test", nodes3d)
        f.read(v, "FileIOTestBase/ThirdRankTensorField3d")
        f.close()
        assert len(v) == len(v0)
        for i in xrange(self.n):
            self.failUnless(v[i] == v0[i],
                            "%s != %s @ %i of %i in ThirdRankTensorField3d test" %
                            (str(v[i]), str(v0[i]), i, self.n))
        self.removeFile("TestThirdRankTensorField3d")
        return

    #---------------------------------------------------------------------------
    # IntFieldList1d
    #---------------------------------------------------------------------------
    def testIntFieldList1d(self):
        fl0 = IntFieldList1d()
        fl0.copyFields()
        fl0.appendNewField("int field 1d control", nodes1d, 0)
        for i in xrange(self.n):
            fl0[0][i] = g.randint(self.intmin, self.intmax)
        assert len(fl0) == 1
        assert len(fl0[0]) == self.n
        f = self.constructor("TestIntFieldList1d", Write)
        f.write(fl0, "FileIOTestBase/IntFieldList1d")
        f.close()
        f = self.constructor("TestIntFieldList1d", Read)
        fl = IntFieldList1d()
        fl.copyFields()
        f.read(fl, "FileIOTestBase/IntFieldList1d")
        f.close()
        assert len(fl) == len(fl0)
        assert len(fl[0]) == len(fl0[0])
        for i in xrange(self.n):
            self.failUnless(fl[0][i] == fl0[0][i],
                            "%i != %i @ %i of %i in IntFieldList1d test" %
                            (fl[0][i], fl0[0][i], i, self.n))
        self.removeFile("TestIntFieldList1d")
        return

    #---------------------------------------------------------------------------
    # IntFieldList2d
    #---------------------------------------------------------------------------
    def testIntFieldList2d(self):
        fl0 = IntFieldList2d()
        fl0.copyFields()
        fl0.appendNewField("int field 2d control", nodes2d, 0)
        for i in xrange(self.n):
            fl0[0][i] = g.randint(self.intmin, self.intmax)
        assert len(fl0) == 1
        assert len(fl0[0]) == self.n
        f = self.constructor("TestIntFieldList2d", Write)
        f.write(fl0, "FileIOTestBase/IntFieldList2d")
        f.close()
        f = self.constructor("TestIntFieldList2d", Read)
        fl = IntFieldList2d()
        fl.copyFields()
        f.read(fl, "FileIOTestBase/IntFieldList2d")
        f.close()
        assert len(fl) == len(fl0)
        assert len(fl[0]) == len(fl0[0])
        for i in xrange(self.n):
            self.failUnless(fl[0][i] == fl0[0][i],
                            "%i != %i @ %i of %i in IntFieldList2d test" %
                            (fl[0][i], fl0[0][i], i, self.n))
        self.removeFile("TestIntFieldList2d")
        return

    #---------------------------------------------------------------------------
    # IntFieldList3d
    #---------------------------------------------------------------------------
    def testIntFieldList3d(self):
        fl0 = IntFieldList3d()
        fl0.copyFields()
        fl0.appendNewField("int field 3d control", nodes3d, 0)
        for i in xrange(self.n):
            fl0[0][i] = g.randint(self.intmin, self.intmax)
        assert len(fl0) == 1
        assert len(fl0[0]) == self.n
        f = self.constructor("TestIntFieldList3d", Write)
        f.write(fl0, "FileIOTestBase/IntFieldList3d")
        f.close()
        f = self.constructor("TestIntFieldList3d", Read)
        fl = IntFieldList3d()
        fl.copyFields()
        f.read(fl, "FileIOTestBase/IntFieldList3d")
        f.close()
        assert len(fl) == len(fl0)
        assert len(fl[0]) == len(fl0[0])
        for i in xrange(self.n):
            self.failUnless(fl[0][i] == fl0[0][i],
                            "%i != %i @ %i of %i in IntFieldList3d test" %
                            (fl[0][i], fl0[0][i], i, self.n))
        self.removeFile("TestIntFieldList3d")
        return

    #---------------------------------------------------------------------------
    # ScalarFieldList1d
    #---------------------------------------------------------------------------
    def testScalarFieldList1d(self):
        fl0 = ScalarFieldList1d()
        fl0.copyFields()
        fl0.appendNewField("scalar field 1d control", nodes1d, 0.0)
        for i in xrange(self.n):
            fl0[0][i] = g.uniform(self.doublemin, self.doublemax)
        assert len(fl0) == 1
        assert len(fl0[0]) == self.n
        f = self.constructor("TestScalarFieldList1d", Write)
        f.write(fl0, "FileIOTestBase/ScalarFieldList1d")
        f.close()
        f = self.constructor("TestScalarFieldList1d", Read)
        fl = ScalarFieldList1d()
        fl.copyFields()
        f.read(fl, "FileIOTestBase/ScalarFieldList1d")
        f.close()
        assert len(fl) == len(fl0)
        assert len(fl[0]) == len(fl0[0])
        for i in xrange(self.n):
            self.failUnless(fl[0][i] == fl0[0][i],
                            "%g != %g @ %i of %i in ScalarFieldList1d test" %
                            (fl[0][i], fl0[0][i], i, self.n))
        self.removeFile("TestScalarFieldList1d")
        return

    #---------------------------------------------------------------------------
    # ScalarFieldList2d
    #---------------------------------------------------------------------------
    def testScalarFieldList2d(self):
        fl0 = ScalarFieldList2d()
        fl0.copyFields()
        fl0.appendNewField("scalar field 2d control", nodes2d, 0.0)
        for i in xrange(self.n):
            fl0[0][i] = g.uniform(self.doublemin, self.doublemax)
        assert len(fl0) == 1
        assert len(fl0[0]) == self.n
        f = self.constructor("TestScalarFieldList2d", Write)
        f.write(fl0, "FileIOTestBase/ScalarFieldList2d")
        f.close()
        f = self.constructor("TestScalarFieldList2d", Read)
        fl = ScalarFieldList2d()
        fl.copyFields()
        f.read(fl, "FileIOTestBase/ScalarFieldList2d")
        f.close()
        assert len(fl) == len(fl0)
        assert len(fl[0]) == len(fl0[0])
        for i in xrange(self.n):
            self.failUnless(fl[0][i] == fl0[0][i],
                            "%g != %g @ %i of %i in ScalarFieldList2d test" %
                            (fl[0][i], fl0[0][i], i, self.n))
        self.removeFile("TestScalarFieldList2d")
        return

    #---------------------------------------------------------------------------
    # ScalarFieldList3d
    #---------------------------------------------------------------------------
    def testScalarFieldList3d(self):
        fl0 = ScalarFieldList3d()
        fl0.copyFields()
        fl0.appendNewField("scalar field 3d control", nodes3d, 0.0)
        for i in xrange(self.n):
            fl0[0][i] = g.uniform(self.doublemin, self.doublemax)
        assert len(fl0) == 1
        assert len(fl0[0]) == self.n
        f = self.constructor("TestScalarFieldList3d", Write)
        f.write(fl0, "FileIOTestBase/ScalarFieldList3d")
        f.close()
        f = self.constructor("TestScalarFieldList3d", Read)
        fl = ScalarFieldList3d()
        fl.copyFields()
        f.read(fl, "FileIOTestBase/ScalarFieldList3d")
        f.close()
        assert len(fl) == len(fl0)
        assert len(fl[0]) == len(fl0[0])
        for i in xrange(self.n):
            self.failUnless(fl[0][i] == fl0[0][i],
                            "%g != %g @ %i of %i in ScalarFieldList3d test" %
                            (fl[0][i], fl0[0][i], i, self.n))
        self.removeFile("TestScalarFieldList3d")
        return

    #---------------------------------------------------------------------------
    # VectorFieldList1d
    #---------------------------------------------------------------------------
    def testVectorFieldList1d(self):
        fl0 = VectorFieldList1d()
        fl0.copyFields()
        fl0.appendNewField("vector field 1d control", nodes1d, Vector1d.zero)
        for i in xrange(self.n):
            fl0[0][i] = Vector1d(g.uniform(self.doublemin, self.doublemax))
        assert len(fl0) == 1
        assert len(fl0[0]) == self.n
        f = self.constructor("TestVectorFieldList1d", Write)
        f.write(fl0, "FileIOTestBase/VectorFieldList1d")
        f.close()
        f = self.constructor("TestVectorFieldList1d", Read)
        fl = VectorFieldList1d()
        fl.copyFields()
        f.read(fl, "FileIOTestBase/VectorFieldList1d")
        f.close()
        assert len(fl) == len(fl0)
        assert len(fl0[0]) == len(fl0[0])
        for i in xrange(self.n):
            self.failUnless(fl[0][i] == fl0[0][i],
                            "%s != %s @ %i of %i in VectorFieldList1d test" %
                            (str(fl[0][i]), str(fl0[0][i]), i, self.n))
        self.removeFile("TestVectorFieldList1d")
        return

    #---------------------------------------------------------------------------
    # VectorFieldList2d
    #---------------------------------------------------------------------------
    def testVectorFieldList2d(self):
        fl0 = VectorFieldList2d()
        fl0.copyFields()
        fl0.appendNewField("vector field 2d control", nodes2d, Vector2d.zero)
        for i in xrange(self.n):
            fl0[0][i] = Vector2d(g.uniform(self.doublemin, self.doublemax),
                                 g.uniform(self.doublemin, self.doublemax))
        assert len(fl0) == 1
        assert len(fl0[0]) == self.n
        f = self.constructor("TestVectorFieldList2d", Write)
        f.write(fl0, "FileIOTestBase/VectorFieldList2d")
        f.close()
        f = self.constructor("TestVectorFieldList2d", Read)
        fl = VectorFieldList2d()
        fl.copyFields()
        f.read(fl, "FileIOTestBase/VectorFieldList2d")
        f.close()
        assert len(fl) == len(fl0)
        assert len(fl[0]) == len(fl0[0])
        for i in xrange(self.n):
            self.failUnless(fl[0][i] == fl0[0][i],
                            "%s != %s @ %i of %i in VectorFieldList2d test" %
                            (str(fl[0][i]), str(fl0[0][i]), i, self.n))
        self.removeFile("TestVectorFieldList2d")
        return

    #---------------------------------------------------------------------------
    # VectorFieldList3d
    #---------------------------------------------------------------------------
    def testVectorFieldList3d(self):
        fl0 = VectorFieldList3d()
        fl0.copyFields()
        fl0.appendNewField("vector field 3d control", nodes3d, Vector3d.zero)
        for i in xrange(self.n):
            fl0[0][i] = Vector3d(g.uniform(self.doublemin, self.doublemax),
                                 g.uniform(self.doublemin, self.doublemax),
                                 g.uniform(self.doublemin, self.doublemax))
        assert len(fl0) == 1
        assert len(fl0[0]) == self.n
        f = self.constructor("TestVectorFieldList3d", Write)
        f.write(fl0, "FileIOTestBase/VectorFieldList3d")
        f.close()
        f = self.constructor("TestVectorFieldList3d", Read)
        fl = VectorFieldList3d()
        fl.copyFields()
        f.read(fl, "FileIOTestBase/VectorFieldList3d")
        f.close()
        assert len(fl) == len(fl0)
        assert len(fl[0]) == len(fl0[0])
        for i in xrange(self.n):
            self.failUnless(fl[0][i] == fl0[0][i],
                            "%s != %s @ %i of %i in VectorFieldList3d test" %
                            (str(fl[0][i]), str(fl0[0][i]), i, self.n))
        self.removeFile("TestVectorFieldList3d")
        return

    #---------------------------------------------------------------------------
    # TensorFieldList1d
    #---------------------------------------------------------------------------
    def testTensorFieldList1d(self):
        fl0 = TensorFieldList1d()
        fl0.copyFields()
        fl0.appendNewField("vector field 1d control", nodes1d, Tensor1d.zero)
        for i in xrange(self.n):
            fl0[0][i] = Tensor1d(g.uniform(self.doublemin, self.doublemax))
        assert len(fl0) == 1
        assert len(fl0[0]) == self.n
        f = self.constructor("TestTensorFieldList1d", Write)
        f.write(fl0, "FileIOTestBase/TensorFieldList1d")
        f.close()
        f = self.constructor("TestTensorFieldList1d", Read)
        fl = TensorFieldList1d()
        fl.copyFields()
        f.read(fl, "FileIOTestBase/TensorFieldList1d")
        f.close()
        assert len(fl) == len(fl0)
        assert len(fl0[0]) == len(fl0[0])
        for i in xrange(self.n):
            self.failUnless(fl[0][i] == fl0[0][i],
                            "%s != %s @ %i of %i in TensorFieldList1d test" %
                            (str(fl[0][i]), str(fl0[0][i]), i, self.n))
        self.removeFile("TestTensorFieldList1d")
        return

    #---------------------------------------------------------------------------
    # TensorFieldList2d
    #---------------------------------------------------------------------------
    def testTensorFieldList2d(self):
        fl0 = TensorFieldList2d()
        fl0.copyFields()
        fl0.appendNewField("vector field 2d control", nodes2d, Tensor2d.zero)
        for i in xrange(self.n):
            fl0[0][i] = Tensor2d(g.uniform(self.doublemin, self.doublemax),
                                 g.uniform(self.doublemin, self.doublemax),
                                 g.uniform(self.doublemin, self.doublemax),
                                 g.uniform(self.doublemin, self.doublemax))
        assert len(fl0) == 1
        assert len(fl0[0]) == self.n
        f = self.constructor("TestTensorFieldList2d", Write)
        f.write(fl0, "FileIOTestBase/TensorFieldList2d")
        f.close()
        f = self.constructor("TestTensorFieldList2d", Read)
        fl = TensorFieldList2d()
        fl.copyFields()
        f.read(fl, "FileIOTestBase/TensorFieldList2d")
        f.close()
        assert len(fl) == len(fl0)
        assert len(fl[0]) == len(fl0[0])
        for i in xrange(self.n):
            self.failUnless(fl[0][i] == fl0[0][i],
                            "%s != %s @ %i of %i in TensorFieldList2d test" %
                            (str(fl[0][i]), str(fl0[0][i]), i, self.n))
        self.removeFile("TestTensorFieldList2d")
        return

    #---------------------------------------------------------------------------
    # TensorFieldList3d
    #---------------------------------------------------------------------------
    def testTensorFieldList3d(self):
        fl0 = TensorFieldList3d()
        fl0.copyFields()
        fl0.appendNewField("vector field 3d control", nodes3d, Tensor3d.zero)
        for i in xrange(self.n):
            fl0[0][i] = Tensor3d(g.uniform(self.doublemin, self.doublemax),
                                 g.uniform(self.doublemin, self.doublemax),
                                 g.uniform(self.doublemin, self.doublemax),
                                 g.uniform(self.doublemin, self.doublemax),
                                 g.uniform(self.doublemin, self.doublemax),
                                 g.uniform(self.doublemin, self.doublemax),
                                 g.uniform(self.doublemin, self.doublemax),
                                 g.uniform(self.doublemin, self.doublemax),
                                 g.uniform(self.doublemin, self.doublemax))
        assert len(fl0) == 1
        assert len(fl0[0]) == self.n
        f = self.constructor("TestTensorFieldList3d", Write)
        f.write(fl0, "FileIOTestBase/TensorFieldList3d")
        f.close()
        f = self.constructor("TestTensorFieldList3d", Read)
        fl = TensorFieldList3d()
        fl.copyFields()
        f.read(fl, "FileIOTestBase/TensorFieldList3d")
        f.close()
        assert len(fl) == len(fl0)
        assert len(fl[0]) == len(fl0[0])
        for i in xrange(self.n):
            self.failUnless(fl[0][i] == fl0[0][i],
                            "%s != %s @ %i of %i in TensorFieldList3d test" %
                            (str(fl[0][i]), str(fl0[0][i]), i, self.n))
        self.removeFile("TestTensorFieldList3d")
        return

    #---------------------------------------------------------------------------
    # SymTensorFieldList1d
    #---------------------------------------------------------------------------
    def testSymTensorFieldList1d(self):
        fl0 = SymTensorFieldList1d()
        fl0.copyFields()
        fl0.appendNewField("vector field 1d control", nodes1d, SymTensor1d.zero)
        for i in xrange(self.n):
            fl0[0][i] = SymTensor1d(g.uniform(self.doublemin, self.doublemax))
        assert len(fl0) == 1
        assert len(fl0[0]) == self.n
        f = self.constructor("TestSymTensorFieldList1d", Write)
        f.write(fl0, "FileIOTestBase/SymTensorFieldList1d")
        f.close()
        f = self.constructor("TestSymTensorFieldList1d", Read)
        fl = SymTensorFieldList1d()
        fl.copyFields()
        f.read(fl, "FileIOTestBase/SymTensorFieldList1d")
        f.close()
        assert len(fl) == len(fl0)
        assert len(fl0[0]) == len(fl0[0])
        for i in xrange(self.n):
            self.failUnless(fl[0][i] == fl0[0][i],
                            "%s != %s @ %i of %i in SymTensorFieldList1d test" %
                            (str(fl[0][i]), str(fl0[0][i]), i, self.n))
        self.removeFile("TestSymTensorFieldList1d")
        return

    #---------------------------------------------------------------------------
    # SymTensorFieldList2d
    #---------------------------------------------------------------------------
    def testSymTensorFieldList2d(self):
        fl0 = SymTensorFieldList2d()
        fl0.copyFields()
        fl0.appendNewField("vector field 2d control", nodes2d, SymTensor2d.zero)
        for i in xrange(self.n):
            xx = g.uniform(self.doublemin, self.doublemax)
            xy = g.uniform(self.doublemin, self.doublemax)
            yy = g.uniform(self.doublemin, self.doublemax)
            fl0[0][i] = SymTensor2d(xx, xy,
                                    xy, yy)
        assert len(fl0) == 1
        assert len(fl0[0]) == self.n
        f = self.constructor("TestSymTensorFieldList2d", Write)
        f.write(fl0, "FileIOTestBase/SymTensorFieldList2d")
        f.close()
        f = self.constructor("TestSymTensorFieldList2d", Read)
        fl = SymTensorFieldList2d()
        fl.copyFields()
        f.read(fl, "FileIOTestBase/SymTensorFieldList2d")
        f.close()
        assert len(fl) == len(fl0)
        assert len(fl[0]) == len(fl0[0])
        for i in xrange(self.n):
            self.failUnless(fl[0][i] == fl0[0][i],
                            "%s != %s @ %i of %i in SymTensorFieldList2d test" %
                            (str(fl[0][i]), str(fl0[0][i]), i, self.n))
        self.removeFile("TestSymTensorFieldList2d")
        return

    #---------------------------------------------------------------------------
    # SymTensorFieldList3d
    #---------------------------------------------------------------------------
    def testSymTensorFieldList3d(self):
        fl0 = SymTensorFieldList3d()
        fl0.copyFields()
        fl0.appendNewField("vector field 3d control", nodes3d, SymTensor3d.zero)
        for i in xrange(self.n):
            xx = g.uniform(self.doublemin, self.doublemax)
            xy = g.uniform(self.doublemin, self.doublemax)
            xz = g.uniform(self.doublemin, self.doublemax)
            yy = g.uniform(self.doublemin, self.doublemax)
            yz = g.uniform(self.doublemin, self.doublemax)
            zz = g.uniform(self.doublemin, self.doublemax)
            fl0[0][i] = SymTensor3d(xx, xy, xz,
                                    xy, yy, yz,
                                    xz, yz, zz)
        assert len(fl0) == 1
        assert len(fl0[0]) == self.n
        f = self.constructor("TestSymTensorFieldList3d", Write)
        f.write(fl0, "FileIOTestBase/SymTensorFieldList3d")
        f.close()
        f = self.constructor("TestSymTensorFieldList3d", Read)
        fl = SymTensorFieldList3d()
        fl.copyFields()
        f.read(fl, "FileIOTestBase/SymTensorFieldList3d")
        f.close()
        assert len(fl) == len(fl0)
        assert len(fl[0]) == len(fl0[0])
        for i in xrange(self.n):
            self.failUnless(fl[0][i] == fl0[0][i],
                            "%s != %s @ %i of %i in SymTensorFieldList3d test" %
                            (str(fl[0][i]), str(fl0[0][i]), i, self.n))
        self.removeFile("TestSymTensorFieldList3d")
        return

    #---------------------------------------------------------------------------
    # writeObject(int)
    #---------------------------------------------------------------------------
    def testWriteObjectInt(self):
        x0 = g.randint(self.intmin, self.intmax)
        f = self.constructor("TestInt", Write)
        f.writeObject(x0, "FileIOTestBase/TestInt")
        f.close()
        f = self.constructor("TestInt", Read)
        x1 = f.readObject("FileIOTestBase/TestInt")
        f.close()
        self.failUnless(x1 == x0,
                        "%i != %i in int OBJECT test" % (x1, x0))
        self.removeFile("TestInt")
        return

    #---------------------------------------------------------------------------
    # writeObject(string)
    #---------------------------------------------------------------------------
    def testWriteObjectString(self):
        x0 = "abcdefg"
        f = self.constructor("TestObject", Write)
        f.writeObject(x0, "FileIOTestBase/TestObject")
        f.close()
        f = self.constructor("TestObject", Read)
        x1 = f.readObject("FileIOTestBase/TestObject")
        f.close()
        self.failUnless(x1 == x0,
                        "%s != %s in string OBJECT test" % (x1, x0))
        self.removeFile("TestObject")
        return

    #---------------------------------------------------------------------------
    # writeObject(list)
    #---------------------------------------------------------------------------
    def testWriteObjectList(self):
        x0 = [49, 492, 59392, 784761, "ackthpt"]
        f = self.constructor("TestObject", Write)
        f.writeObject(x0, "FileIOTestBase/TestObject")
        f.close()
        f = self.constructor("TestObject", Read)
        x1 = f.readObject("FileIOTestBase/TestObject")
        f.close()
        self.failUnless(x1 == x0,
                        "%s != %s in list OBJECT test" % (x1, x0))
        self.removeFile("TestObject")
        return

    #---------------------------------------------------------------------------
    # writeObject(Vector3d)
    #---------------------------------------------------------------------------
    def testWriteObjectVector3d(self):
        x0 = Vector3d(g.uniform(self.doublemin, self.doublemax),
                      g.uniform(self.doublemin, self.doublemax),
                      g.uniform(self.doublemin, self.doublemax))
        x1 = Vector3d()
        f = self.constructor("TestVector3d", Write)
        f.writeObject(x0, "FileIOTestBase/TestVector3d")
        f.close()
        f = self.constructor("TestVector3d", Read)
        x1 = f.readObject("FileIOTestBase/TestVector3d")
        f.close()
        self.failUnless(x1 == x0,
                        "%s != %s in Vector3d OBJECT test" % (str(x1), str(x0)))
        self.removeFile("TestVector3d")
        return
