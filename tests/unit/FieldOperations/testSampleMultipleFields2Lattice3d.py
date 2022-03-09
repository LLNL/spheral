#ATS:test(SELF,       label="sampleMultipleFields2Lattice 3D unit tests (serial)")
#ATS:test(SELF, np=4, label="sampleMultipleFields2Lattice 3D unit tests (4 proc)")

import unittest

from Spheral import *
from SpheralTestUtilities import fuzzyEqual
import mpi

from testSampleMultipleFields2Lattice import TestSampleMultipleFields2Lattice

#===============================================================================
# 2-D tests.
#===============================================================================
class TestSampleMultipleFields2Lattice3d(TestSampleMultipleFields2Lattice,
                                         unittest.TestCase):

    #---------------------------------------------------------------------------
    # Initialize the problem.
    #---------------------------------------------------------------------------
    def setUp(self):

        from Spheral3d import (vector_of_int, Vector, Tensor, GammaLawGasMKS,
                               TableKernel, BSplineKernel, makeFluidNodeList,
                               ScalarField, VectorField, DataBase, Plane,
                               PeriodicBoundary)

        self.ndim = 3
        self.genxmin = (0.0, 0.0, 0.0)
        self.genxmax = (1.0, 1.0, 1.0)
        self.xmin = Vector(0.2, 0.2, 0.2)
        self.xmax = Vector(0.8, 0.8, 0.8)
        self.nsample = vector_of_int()
        [self.nsample.append(x) for x in (40, 40, 40)]

        # Tolerances for the test
        self.scalarTol = 1.0e-2
        self.vectorTol = 1.0e-2
        self.tensorTol = 1.0e-2

        nx, ny, nz = 20, 20, 20
        self.rho0 = 10.0
        self.v0 = Vector(1.0, -1.0, 0.5)
        self.eps0 = -1.0
        self.gradv0 = Tensor(8.5, -4.0, 10.0,
                             2.2,  1.3, -1.0,
                             5.1, -2.4, 14.5)

        # Create the nodes and such.
        self.eos = GammaLawGasMKS(5.0/3.0, 1.0)
        self.WT = TableKernel(BSplineKernel())
        self.nodes = makeFluidNodeList("nodes", self.eos)

        # Distribute the nodes.
        from GenerateNodeDistribution3d import GenerateNodeDistribution3d as GenerateNodeDistribution
        from PeanoHilbertDistributeNodes import distributeNodes3d as distributeNodes
        generator = GenerateNodeDistribution(nx, ny, nz,
                                             self.rho0,
                                             "lattice",
                                             xmin = self.genxmin,
                                             xmax = self.genxmax,
                                             nNodePerh = 2.01)
        distributeNodes((self.nodes, generator))

        # Set the velocities and energies.
        self.nodes.velocity(VectorField("tmp", self.nodes, self.v0))
        self.nodes.specificThermalEnergy(ScalarField("tmp", self.nodes, self.eps0))

        self.db = DataBase()
        self.db.appendNodeList(self.nodes)

        # Create the boundary conditions.
        px0 = Plane(Vector(0, 0, 0), Vector( 1,  0,  0))
        px1 = Plane(Vector(1, 0, 0), Vector(-1,  0,  0))
        py0 = Plane(Vector(0, 0, 0), Vector( 0,  1,  0))
        py1 = Plane(Vector(0, 1, 0), Vector( 0, -1,  0))
        pz0 = Plane(Vector(0, 0, 0), Vector( 0,  0,  1))
        pz1 = Plane(Vector(0, 0, 1), Vector( 0,  0, -1))
        xbc = PeriodicBoundary(px0, px1)
        ybc = PeriodicBoundary(py0, py1)
        zbc = PeriodicBoundary(pz0, pz1)
        self.bcs = [xbc, ybc, zbc]
        try:
            self.bcs.append(TreeDistributedBoundary3d.instance())
        except:
            if mpi.procs > 1:
                raise RuntimeError, "Unable to get parallel boundary condition"
            else:
                pass

        # Enforce boundaries.
        db = DataBase()
        db.appendNodeList(self.nodes)
        for bc in self.bcs:
            bc.setAllGhostNodes(db)
            bc.finalizeGhostBoundary()
            self.nodes.neighbor().updateNodes()
        for bc in self.bcs:
            bc.applyGhostBoundary(self.nodes.mass())
            bc.applyGhostBoundary(self.nodes.massDensity())
            bc.applyGhostBoundary(self.nodes.specificThermalEnergy())
            bc.applyGhostBoundary(self.nodes.velocity())
        for bc in self.bcs:
            bc.finalizeGhostBoundary()

        self.H0 = self.nodes.Hfield()[0]
        return

    def tearDown(self):
        del self.nodes

#===============================================================================
# Run the tests
#===============================================================================
if __name__ == "__main__":
    unittest.main()
