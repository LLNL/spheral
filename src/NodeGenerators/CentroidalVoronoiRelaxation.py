# Apply various centroidal Voronoi relaxation methods to optimize NodeGenerators.

from Spheral2d import *
from generateMesh import *

#-------------------------------------------------------------------------------
# Relax nodes on fixed radii.
#-------------------------------------------------------------------------------
class RadialCentroidalRelaxation:
    def __init__(self, origin,
                 tolerance = 1.0e-5,
                 maxIterations = 100):
        self.origin = origin
        self.tolerance = tolerance
        self.maxIterations = maxIterations
        return

    def __call__(self, generator):
        n = generator.localNumNodes()
        nodes = makeVoidNodeList("temp void nodes",
                                 numInternal = n)
        pos = nodes.positions()
        rnodes = []
        for i in range(n):
            pos[i] = Vector(generator.x[i], generator.y[i])
            rnodes.append((pos[i] - self.origin).magnitude())
        assert len(rnodes) == n
        if generator.xmin:
            xmin = Vector(generator.xmin[0], generator.xmin[1])
        else:
            xmin = None
        if generator.xmax:
            xmax = Vector(generator.xmax[0], generator.xmax[1])
        else:
            xmax = None

        # Iterate until we either converge or hit the max iterations.
        maxDelta = 10.0*self.tolerance
        iter = 0
        while maxDelta > self.tolerance and iter < self.maxIterations:
            maxDelta = 0.0
            iter += 1
            mesh, void = generatePolygonalMesh([nodes], 
                                               xmin = None, # xmin,
                                               xmax = None, # xmax,
                                               generateVoid = False,
                                               removeBoundaryZones = True)
            assert mesh.numZones == n
            for izone in range(n):
                zone = mesh.zone(izone)
                centroid = zone.position()
                rhat = (centroid - self.origin).unitVector()
                newpos = rnodes[izone]*rhat
                maxDelta = max(maxDelta, (newpos - pos[izone]).magnitude())
                pos[izone] = newpos
            maxDelta = mpi.allreduce(maxDelta, mpi.MAX)
            print("Iteration %i : max delta = %g" % (iter, maxDelta))

        # Assign the postions.
        generator.x = [pos[i].x for i in range(nodes.numInternalNodes)]
        generator.y = [pos[i].y for i in range(nodes.numInternalNodes)]

        return
