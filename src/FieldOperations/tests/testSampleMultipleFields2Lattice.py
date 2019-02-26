#ATS:test(SELF, label="sampleMultipleFields2Lattice unit tests (serial)")
#ATS:test(SELF, np=4, label="sampleMultipleFields2Lattice unit tests (parallel)")

import unittest

from Spheral import *
from SpheralTestUtilities import fuzzyEqual
import mpi

#===============================================================================
# Base test class -- generic for all dimensions.
#===============================================================================
class TestSampleMultipleFields2Lattice:

    #---------------------------------------------------------------------------
    # Run the test!
    #---------------------------------------------------------------------------
    def testSample(self):

        if self.ndim == 1:
            from Spheral1d import (ScalarFieldList, VectorFieldList, TensorFieldList, SymTensorFieldList, 
                                   FieldListSet, sampleMultipleFields2LatticeMash,
                                   vector_of_vector_of_Vector, vector_of_vector_of_Tensor, vector_of_vector_of_SymTensor)
        elif self.ndim == 2:
            from Spheral2d import (ScalarFieldList, VectorFieldList, TensorFieldList, SymTensorFieldList, 
                                   FieldListSet, sampleMultipleFields2LatticeMash,
                                   vector_of_vector_of_Vector, vector_of_vector_of_Tensor, vector_of_vector_of_SymTensor)
        else:
            assert ndim == 3
            from Spheral3d import (ScalarFieldList, VectorFieldList, TensorFieldList, SymTensorFieldList, 
                                   FieldListSet, sampleMultipleFields2LatticeMash,
                                   vector_of_vector_of_Vector, vector_of_vector_of_Tensor, vector_of_vector_of_SymTensor)

        # Build the FieldLists we're going to sample.
        rho = ScalarFieldList()
        eps = ScalarFieldList()
        vel = VectorFieldList()
        Hfl = SymTensorFieldList()
        rho.appendField(self.nodes.massDensity())
        eps.appendField(self.nodes.specificThermalEnergy())
        vel.appendField(self.nodes.velocity())
        Hfl.appendField(self.nodes.Hfield())

        # Put them into a FieldListSet.
        fieldListSet = FieldListSet()
        fieldListSet.ScalarFieldLists.append(rho)
        fieldListSet.ScalarFieldLists.append(eps)
        fieldListSet.VectorFieldLists.append(vel)
        fieldListSet.SymTensorFieldLists.append(Hfl)

        # Build the mask.
        mask = self.db.newGlobalIntFieldList(1)

        # Now sample those suckers.
        r = VectorFieldList()
        w = ScalarFieldList()
        H = SymTensorFieldList()
        r.appendField(self.nodes.positions())
        w.appendField(self.nodes.mass())
        H.appendField(self.nodes.Hfield())
        scalar_samples, vector_samples, tensor_samples, symtensor_samples = sampleMultipleFields2LatticeMash(fieldListSet,
                                                                                                             r,
                                                                                                             w,
                                                                                                             H,
                                                                                                             mask,
                                                                                                             self.WT,
                                                                                                             self.xmin,
                                                                                                             self.xmax,
                                                                                                             self.nsample)

        # Did we get back the correct numbers of sampled values?
        assert len(scalar_samples) == 2
        assert len(vector_samples) == 1
        assert len(tensor_samples) == 0
        assert len(symtensor_samples) == 1
        ntot = 1
        for x in self.nsample:
            ntot *= x
        for array in scalar_samples:
            assert mpi.allreduce(len(array), mpi.SUM) == ntot
        for array in vector_samples:
            assert mpi.allreduce(len(array), mpi.SUM) == ntot
        for array in tensor_samples:
            assert mpi.allreduce(len(array), mpi.SUM) == ntot
        for array in symtensor_samples:
            assert mpi.allreduce(len(array), mpi.SUM) == ntot
        
        # Extract the arrays.
        rhoarray = scalar_samples[0]
        epsarray = scalar_samples[1]
        velarray = vector_samples[0]
        Harray = symtensor_samples[0]

        # See if the resulting fields are constant enough.
        self.failUnless(min([fuzzyEqual(x, self.rho0, self.scalarTol) for x in rhoarray]),
                        "Failing rho comparison: expect=%g min=%g, max=%g, tol=%g" % (self.rho0,
                                                                                      min(rhoarray),
                                                                                      max(rhoarray),
                                                                                      self.scalarTol))
        self.failUnless(min([fuzzyEqual(x, self.eps0, self.scalarTol) for x in epsarray]),
                        "Failing eps comparison: expect=%g min=%g, max=%g, tol=%g" % (self.eps0,
                                                                                      min(epsarray),
                                                                                      max(epsarray),
                                                                                      self.scalarTol))
        self.failUnless(min([fuzzyEqual((x - self.v0).magnitude(), 0.0, self.vectorTol) for x in velarray]),
                        "Failing vel comparison: min=%g, max=%g, tol=%g" % (min((x - self.v0).magnitude() for x in velarray),
                                                                            max((x - self.v0).magnitude() for x in velarray),
                                                                            self.vectorTol))
##         self.failUnless(min([fuzzyEqual((x - self.gradv0).doubledot(x - self.gradv0), 0.0, self.tensorTol) for x in gradvelarray]),
##                         "Failing gradv comparison: min=%g, max=%g, tol=%g" % (min([(x - self.gradv0).doubledot(x - self.gradv0) for x in gradvelarray]),
##                                                                               max([(x - self.gradv0).doubledot(x - self.gradv0) for x in gradvelarray]),
##                                                                               self.tensorTol))
        self.failUnless(min([fuzzyEqual((x - self.H0).doubledot(x - self.H0), 0.0, self.tensorTol) for x in Harray]),
                        "Failing H comparison: min=%g, max=%g, tol=%g" % (min([(x - self.H0).doubledot(x - self.H0) for x in Harray]),
                                                                          max([(x - self.H0).doubledot(x - self.H0) for x in Harray]),
                                                                          self.tensorTol))

        return

#===============================================================================
# 1-D tests.
#===============================================================================
class TestSampleMultipleFields2Lattice1d(TestSampleMultipleFields2Lattice,
                                         unittest.TestCase):

    #---------------------------------------------------------------------------
    # Initialize the problem.
    #---------------------------------------------------------------------------
    def setUp(self):

        self.ndim = 1
        self.xmin = Vector1d(0.0)
        self.xmax = Vector1d(1.0)
        self.nsample = vector_of_int()
        self.nsample.append(100)

        # Tolerances for the test
        self.scalarTol = 1.0e-5
        self.vectorTol = 1.0e-3
        self.tensorTol = 1.0e-4

        n = 100
        self.rho0 = 10.0
        self.v0 = Vector1d(1.0)
        self.eps0 = -1.0
        self.gradv0 = Tensor1d(8.0)
        x0, x1 = 0.0, 1.0
        
        # Create the nodes and such.
        self.eos = GammaLawGasMKS1d(5.0/3.0, 1.0)
        self.WT = TableKernel1d(BSplineKernel1d())
        self.nodes = makeFluidNodeList1d("nodes", self.eos)
        self.neighbor = self.nodes.neighbor()

        # Distribute the nodes.
        from DistributeNodes import distributeNodesInRange1d
        distributeNodesInRange1d([(self.nodes, n, self.rho0, (x0, x1))])

        # Set the velocities and energies.
        self.nodes.velocity(VectorField1d("tmp", self.nodes, self.v0))
        self.nodes.specificThermalEnergy(ScalarField1d("tmp", self.nodes, self.eps0))

        self.db = DataBase1d()
        self.db.appendNodeList(self.nodes)

        # Create the boundary conditions.
        p0 = Plane1d(Vector1d(0.0), Vector1d(1.0))
        p1 = Plane1d(Vector1d(1.0), Vector1d(-1.0))
        xbc = PeriodicBoundary1d(p0, p1)
        self.bcs = [xbc]
        try:
            dbc = TreeDistributedBoundary1d.instance()
            self.bcs.append(dbc)
        except:
            if mpi.procs > 1:
                raise RuntimeError, "Unable to get parallel boundary condition"
            else:
                pass

        # Enforce boundaries.
        db = DataBase1d()
        db.appendNodeList(self.nodes)
        for bc in self.bcs:
            bc.setAllGhostNodes(db)
            bc.finalizeGhostBoundary()
            self.neighbor.updateNodes()
        for bc in self.bcs:
            bc.applyGhostBoundary(self.nodes.mass())
            bc.applyGhostBoundary(self.nodes.massDensity())
            bc.applyGhostBoundary(self.nodes.specificThermalEnergy())
            bc.applyGhostBoundary(self.nodes.velocity())
        for bc in self.bcs:
            bc.finalizeGhostBoundary()

        self.H0 = self.nodes.Hfield()[0]
        return

#===============================================================================
# 2-D tests.
#===============================================================================
class TestSampleMultipleFields2Lattice2d(TestSampleMultipleFields2Lattice,
                                         unittest.TestCase):

    #---------------------------------------------------------------------------
    # Initialize the problem.
    #---------------------------------------------------------------------------
    def setUp(self):

        self.ndim = 2
        self.genxmin = (0.0, 0.0)
        self.genxmax = (1.0, 1.0)
        self.xmin = Vector2d(0.2, 0.2)
        self.xmax = Vector2d(0.8, 0.8)
        self.nsample = vector_of_int()
        [self.nsample.append(x) for x in (100, 100)]

        # Tolerances for the test
        self.scalarTol = 1.0e-2
        self.vectorTol = 1.0e-2
        self.tensorTol = 1.0e-2

        nx, ny = 50, 50
        self.rho0 = 10.0
        self.v0 = Vector2d(1.0, -1.0)
        self.eps0 = -1.0
        self.gradv0 = Tensor2d(8.5, -4.0,
                               2.2, 1.3)

        # Create the nodes and such.
        self.eos = GammaLawGasMKS2d(5.0/3.0, 1.0)
        self.WT = TableKernel2d(BSplineKernel2d())
        self.nodes = makeFluidNodeList2d("nodes", self.eos)
        self.neighbor = self.nodes.neighbor()

        # Distribute the nodes.
        from GenerateNodeDistribution2d import GenerateNodeDistribution2d
        from DistributeNodes import distributeNodes2d
        generator = GenerateNodeDistribution2d(nx, ny,
                                               self.rho0,
                                               "lattice",
                                               xmin = self.genxmin,
                                               xmax = self.genxmax,
                                               nNodePerh = 2.01)
        distributeNodes2d((self.nodes, generator))

        # Set the velocities and energies.
        self.nodes.velocity(VectorField2d("tmp", self.nodes, self.v0))
        self.nodes.specificThermalEnergy(ScalarField2d("tmp", self.nodes, self.eps0))

        self.db = DataBase2d()
        self.db.appendNodeList(self.nodes)

        # Create the boundary conditions.
        px0 = Plane2d(Vector2d(0.0, 0.0), Vector2d(1.0, 0.0))
        px1 = Plane2d(Vector2d(1.0, 0.0), Vector2d(-1.0, 0.0))
        py0 = Plane2d(Vector2d(0.0, 0.0), Vector2d(0.0, 1.0))
        py1 = Plane2d(Vector2d(0.0, 1.0), Vector2d(0.0, -1.0))
        xbc = PeriodicBoundary2d(px0, px1)
        ybc = PeriodicBoundary2d(py0, py1)
        self.bcs = [xbc, ybc]
        try:
            dbc = TreeDistributedBoundary2d.instance()
            self.bcs.append(dbc)
        except:
            if mpi.procs > 1:
                raise RuntimeError, "Unable to get parallel boundary condition"
            else:
                pass

        # Enforce boundaries.
        db = DataBase2d()
        db.appendNodeList(self.nodes)
        for bc in self.bcs:
            bc.setAllGhostNodes(db)
            bc.finalizeGhostBoundary()
            self.neighbor.updateNodes()
        for bc in self.bcs:
            bc.applyGhostBoundary(self.nodes.mass())
            bc.applyGhostBoundary(self.nodes.massDensity())
            bc.applyGhostBoundary(self.nodes.specificThermalEnergy())
            bc.applyGhostBoundary(self.nodes.velocity())
        for bc in self.bcs:
            bc.finalizeGhostBoundary()

        self.H0 = self.nodes.Hfield()[0]
        return

#===============================================================================
# 3-D tests.
#===============================================================================
## class TestSampleMultipleFields2Lattice3d(TestSampleMultipleFields2Lattice,
##                                          unittest.TestCase):

##     #---------------------------------------------------------------------------
##     # Initialize the problem.
##     #---------------------------------------------------------------------------
##     def setUp(self):

##         self.ScalarFieldList = ScalarFieldList3d
##         self.VectorFieldList = VectorFieldList3d
##         self.TensorFieldList = TensorFieldList3d
##         self.SymTensorFieldList = SymTensorFieldList3d
##         self.FieldListSet = FieldListSet3d
##         self.genxmin = (0.0, 0.0, 0.0)
##         self.genxmax = (1.0, 1.0, 1.0)
##         self.xmin = Vector3d(0.2, 0.2, 0.2)
##         self.xmax = Vector3d(0.8, 0.8, 0.8)
##         self.nsample = (100, 100, 100)
##         self.sample = sampleMultipleFields2LatticeMash3d

##         # Tolerances for the test
##         self.scalarTol = 1.0e-3
##         self.vectorTol = 1.0e-3
##         self.tensorTol = 1.0e-3

##         nx, ny, nz = 50, 50, 50
##         self.rho0 = 10.0
##         self.v0 = Vector3d(1.0, -1.0, 0.5)
##         self.eps0 = -1.0
##         self.gradv0 = Tensor3d( 8.5, -4.0, -2.0,
##                                 2.2,  1.3,  3.5,
##                                -1.8,  2.0,  1.0)
        
##         # Create the nodes and such.
##         self.eos = GammaLawGasMKS3d(5.0/3.0, 1.0)
##         self.WT = TableKernel3d(BSplineKernel3d())
##         self.nodes = makeFluidNodeList3d("nodes", self.eos)
##         self.neighbor = self.nodes.neighbor()

##         # Distribute the nodes.
##         from GenerateNodeDistribution3d import GenerateNodeDistribution3d
##         from DistributeNodes import distributeNodes3d
##         generator = GenerateNodeDistribution3d(nx, ny, nz,
##                                                self.rho0,
##                                                "lattice",
##                                                xmin = self.genxmin,
##                                                xmax = self.genxmax,
##                                                nNodePerh = 2.01)
##         distributeNodes3d((self.nodes, generator))
##         self.H0 = self.nodes.Hfield()[0]

##         # Set the velocities and energies.
##         self.nodes.velocity(VectorField3d("tmp", self.nodes, self.v0))
##         self.nodes.specificThermalEnergy(ScalarField3d("tmp", self.nodes, self.eps0))

##         return

#===============================================================================
# Run the tests
#===============================================================================
if __name__ == "__main__":
    unittest.main()
