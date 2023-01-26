from math import *
import mpi

from NodeGeneratorBase import *
from Spheral import Vector2d, Tensor2d, SymTensor2d, \
     rotationMatrix2d, testPointInBox2d

import PolytopeModules as poly
from sobol import i4_sobol

#-------------------------------------------------------------------------------
# 2D Generator.  Seeds positions within the given boundary using the Sobol 
# sequence to create the seeds.
#-------------------------------------------------------------------------------
class CentroidalGenerator2d(NodeGeneratorBase):

    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self,
                 n,
                 rho,
                 boundary,
                 holes = [],
                 maxIterations = 100,
                 fracTol = 1.0e-3,
                 tessellationFileName = None,
                 nNodePerh = 2.01,
                 offset = (0.0, 0.0),
                 rejecter = None):

        assert n > 0

        # Did we get passed a function or a constant for the density?
        if type(rho) == type(1.0):
            def rhofunc(posi):
                return rho
        else:
            rhofunc = rho
        self.rhofunc = rhofunc

        # Build a polytope PLC version of the boundary.
        plc = poly.polytope.PLC2d()
        plc_coords = poly.vector_of_double()
        edges = boundary.edges
        vertices = boundary.vertices()
        plc.facets.resize(len(edges))
        for i, edge in enumerate(edges):
            plc.facets[i].append(edge.first)
            plc.facets[i].append(edge.second)
            assert len(plc.facets[i]) == 2
        for p in vertices:
            plc_coords.append(p[0])
            plc_coords.append(p[1])
        assert len(plc_coords) == 2*len(vertices)

        # Add any holes to the boundary PLC.
        plc.holes.resize(len(holes))
        for ihole, hole in enumerate(holes):
            offlast = len(plc_coords)/2
            edges = hole.edges
            vertices = hole.vertices()
            plc.holes[ihole].resize(len(edges))
            for i, edge in enumerate(edges):
                plc.holes[ihole][i].append(offlast + edge.first)
                plc.holes[ihole][i].append(offlast + edge.second)
                assert len(plc.holes[ihole][i]) == 2
            for p in vertices:
                plc_coords.append(p[0])
                plc_coords.append(p[1])
            assert len(plc_coords) % 2 == 0

        # Initialize the desired number of generators in the boundary using the Sobol sequence.
        generators = poly.vector_of_double()
        seed = 0
        length = max(boundary.xmax.x - boundary.xmin.x,
                     boundary.xmax.y - boundary.xmin.y)
        while len(generators) < 2*n:
            [coords, seed] = i4_sobol(2, seed)
            p = boundary.xmin + length*Vector2d(coords[0], coords[1])
            ihole = 0
            use = boundary.contains(p, False)
            if use:
                while use and ihole < len(holes):
                    use = not holes[ihole].contains(p, False)
                    ihole += 1
            if use:
                generators.append(p.x)
                generators.append(p.y)
        assert len(generators) == 2*n

        # Iterate the points toward centroidal relaxation.
        self.tessellation = poly.polytope.Tessellation2d()
        tessellator = poly.polytope.BoostTessellator2d()
        iteration = 0
        maxDelta = 2.0*fracTol
        while iteration < maxIterations and maxDelta > fracTol:
            tessellator.tessellate(points = generators,
                                   PLCpoints = plc_coords,
                                   geometry = plc,
                                   mesh = self.tessellation)
            new_generators = self.computeWeightedCentroids(self.tessellation)
            assert len(new_generators) == len(generators)
            maxDelta = 0.0
            for i in range(len(generators)/2):
                deltai = sqrt((generators[2*i] - new_generators[2*i])**2 +
                              (generators[2*i+1] - new_generators[2*i+1])**2)
                maxDelta = max(maxDelta, deltai/length)
                generators[2*i] = 0.5*(generators[2*i] + new_generators[2*i])
                generators[2*i+1] = 0.5*(generators[2*i+1] + new_generators[2*i+1])
            iteration += 1
            print("CentroidalGenerator2d: Iteration %i, maxDelta=%g" % (iteration, maxDelta))

        # If requested, write out the final tessellation to a silo file.
        if tessellationFileName:
            poly.polytope.writeTessellation2d(mesh = self.tessellation,
                                              filePrefix = tessellationFileName,
                                              nodeFields = None,
                                              edgeFields = None,
                                              faceFields = None,
                                              cellFields = None,
                                              cycle = iteration)

        # Now we can fill out the usual Spheral generator info.
        assert len(self.tessellation.cells) == n
        self.x, self.y, self.m, self.H = [], [], [], []
        centroids = self.computeWeightedCentroids(self.tessellation)
        masses = self.computeMasses(self.tessellation)
        areas = self.computeAreas(self.tessellation)
        assert len(centroids) == 2*n
        assert len(masses) == n
        assert len(areas) == n
        for i in range(n):
            self.x.append(centroids[2*i] + offset[0])
            self.y.append(centroids[2*i+1] + offset[1])
            self.m.append(masses[i])
            hi = nNodePerh * sqrt(areas[i]/pi)
            assert hi > 0.0
            self.H.append(SymTensor2d(1.0/hi, 0.0, 0.0, 1.0/hi))
        assert len(self.x) == n
        assert len(self.y) == n
        assert len(self.m) == n
        assert len(self.H) == n

        # If the user provided a "rejecter", give it a pass
        # at the nodes.
        if rejecter:
            self.x, self.y, self.m, self.H = rejecter(self.x,
                                                      self.y,
                                                      self.m,
                                                      self.H)

        # Have the base class break up the serial node distribution
        # for parallel cases.
        NodeGeneratorBase.__init__(self, True,
                                   self.x, self.y, self.m, self.H)
        return
            

    #---------------------------------------------------------------------------
    # Compute the barycenters of a tessellation.
    #---------------------------------------------------------------------------
    def computeBarycenters(self, tessellation):
        xb       = [0.0 for i in range(tessellation.cells.size())]
        yb       = [0.0 for i in range(tessellation.cells.size())]
        numfaces = [cell.size() for cell in tessellation.cells]
        for iface,face in enumerate(tessellation.faces):
            fc = tessellation.faceCells[iface]
            xf = 0.5*(tessellation.nodes[2*face[0]  ] + tessellation.nodes[2*face[1]  ])
            yf = 0.5*(tessellation.nodes[2*face[0]+1] + tessellation.nodes[2*face[1]+1])
            for c in fc:
                if c < 0: icell = ~c
                else:     icell = c
                xb[icell] += xf
                yb[icell] += yf
        return [(x/n, y/n) for x,y,n in zip(xb,yb,numfaces)]
       
    #---------------------------------------------------------------------------
    # Compute the areas of a tessellation.
    #---------------------------------------------------------------------------
    def computeAreas(self, tessellation):
        result = []
        barycenters = self.computeBarycenters(tessellation)
        for icell, cell in enumerate(tessellation.cells):
            area = 0.0
            for ftmp in cell:
                if ftmp < 0:  
                    n0 = tessellation.faces[~ftmp][1]
                    n1 = tessellation.faces[~ftmp][0]
                else:         
                    n0 = tessellation.faces[ ftmp][0]
                    n1 = tessellation.faces[ ftmp][1]
                x0 = tessellation.nodes[2*n0  ]
                y0 = tessellation.nodes[2*n0+1]
                x1 = tessellation.nodes[2*n1  ]
                y1 = tessellation.nodes[2*n1+1]
                xb = barycenters[icell][0]
                yb = barycenters[icell][1]
                xt = (x0 + x1 + xb)/3.0
                yt = (y0 + y1 + yb)/3.0
                xe = (x0 + x1)/2.0
                ye = (y0 + y1)/2.0
                area += (x0-xb)*(y1-yb) - (x1-xb)*(y0-yb)
            result.append(0.5*area)
        return result

    #---------------------------------------------------------------------------
    # Compute the masses a tessellation.
    #---------------------------------------------------------------------------
    def computeMasses(self, tessellation):
        barycenters = self.computeBarycenters(tessellation)
        result = []
        for icell, cell in enumerate(tessellation.cells):
            mass = 0.0
            for ftmp in cell:
                if ftmp < 0:  
                    n0 = tessellation.faces[~ftmp][1]
                    n1 = tessellation.faces[~ftmp][0]
                else:         
                    n0 = tessellation.faces[ ftmp][0]
                    n1 = tessellation.faces[ ftmp][1]
                x0 = tessellation.nodes[2*n0  ]
                y0 = tessellation.nodes[2*n0+1]
                x1 = tessellation.nodes[2*n1  ]
                y1 = tessellation.nodes[2*n1+1]
                xb = barycenters[icell][0]
                yb = barycenters[icell][1]
                xt = (x0 + x1 + xb)/3.0
                yt = (y0 + y1 + yb)/3.0
                xe = (x0 + x1)/2.0
                ye = (y0 + y1)/2.0
                d  = 0.5*((x0-xb)*(y1-yb) - (x1-xb)*(y0-yb))*self.rhofunc(Vector2d(xe, ye))
                mass += d
            result.append(mass)
        return result

    #---------------------------------------------------------------------------
    # Compute the weighted centroids of a tessellation.
    #---------------------------------------------------------------------------
    def computeWeightedCentroids(self, tessellation):
        barycenters = self.computeBarycenters(tessellation)
        centroids = []
        for icell, cell in enumerate(tessellation.cells):
            mass = 0.0
            xc = 0.0
            yc = 0.0
            for ftmp in cell:
                if ftmp < 0:  
                    n0 = tessellation.faces[~ftmp][1]
                    n1 = tessellation.faces[~ftmp][0]
                else:         
                    n0 = tessellation.faces[ ftmp][0]
                    n1 = tessellation.faces[ ftmp][1]
                x0 = tessellation.nodes[2*n0  ]
                y0 = tessellation.nodes[2*n0+1]
                x1 = tessellation.nodes[2*n1  ]
                y1 = tessellation.nodes[2*n1+1]
                xb = barycenters[icell][0]
                yb = barycenters[icell][1]
                xt = (x0 + x1 + xb)/3.0
                yt = (y0 + y1 + yb)/3.0
                xe = (x0 + x1)/2.0
                ye = (y0 + y1)/2.0
                d  = 0.5*((x0-xb)*(y1-yb) - (x1-xb)*(y0-yb))*self.rhofunc(Vector2d(xe, ye))
                mass += d
                xc   += d*xt
                yc   += d*yt
            xc   /= mass
            yc   /= mass
            centroids.append(xc)
            centroids.append(yc)
        return centroids

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        assert i >= 0 and i < len(self.x)
        assert len(self.x) == len(self.y)
        return Vector2d(self.x[i], self.y[i])
    
    #---------------------------------------------------------------------------
    # Get the mass for the given node index.
    #---------------------------------------------------------------------------
    def localMass(self, i):
        assert i >= 0 and i < len(self.m)
        return self.m[i]
    
    #---------------------------------------------------------------------------
    # Get the mass density for the given node index.
    #---------------------------------------------------------------------------
    def localMassDensity(self, i):
        assert i >= 0 and i < len(self.H)
        return self.rhofunc(Vector2d(self.x[i], self.y[i]))
    
    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]
