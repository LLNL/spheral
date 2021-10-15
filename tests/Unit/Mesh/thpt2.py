from Spheral2d import *
from siloMeshDump import *
from generateMesh import *
from math import *
import random

eos = GammaLawGasMKS(5.0/3.0, 1.0)
xbc = ReflectingBoundary(Plane(Vector(-100,-100), Vector(1,0)))
ybc = ReflectingBoundary(Plane(Vector(-100,-100), Vector(0,1)))
bcs = vector_of_Boundary()
bcs.append(xbc)
bcs.append(ybc)

W = TableKernel(BSplineKernel(), 1000)

def createNodes(gens, dx):
    nodes = makeFluidNodeList("some_nodes", eos,
                              numInternal = len(gens))
    pos = nodes.positions()
    H = nodes.Hfield()
    mass = nodes.mass()
    rho = nodes.massDensity()
    vel = nodes.velocity()
    H0 = 1.0/(dx*nodes.nodesPerSmoothingScale) * SymTensor.one
    for i in xrange(len(gens)):
        xi = gens[i]
        pos[i] = xi
        H[i] = H0
        mass[i] = 1.0
        rho[i] = 1.0 + xi.magnitude2()
        vel[i] = xi

    db = DataBase()
    db.appendNodeList(nodes)
    for bc in (xbc, ybc):
        bc.setAllGhostNodes(db)
        bc.finalizeGhostBoundary()
        nodes.neighbor().updateNodes()
    db.updateConnectivityMap()

    iterateIdealH(db, bcs, W, SPHSmoothingScale(),
                  tolerance = 1.0e-4)
    db.updateConnectivityMap()

    return nodes, db

generators = vector_of_Vector()
nx = 30
dx = 1.0
print "Creating generators for regular mesh."
for iy in xrange(nx):
    for ix in xrange(nx):
        generators.append(Vector((ix + 0.5)*dx,
                                 (iy + 0.5)*dx))
nodes, db = createNodes(generators, dx)
print "Generating regular mesh."
mesh, void = generatePolygonalMesh([nodes], bcs,
                                   generateVoid = False,
                                   removeBoundaryZones = True)

## mesh = PolygonalMesh(generators, Vector(0,0), Vector(nx,nx))
## print "Writing..."
## siloMeshDump("testPolygonalQuads", mesh, nodeLists = [nodes])
## print "Done."
## del nodes

## generators = vector_of_Vector()
## n = nx**2
## print "Creating generators for random mesh."
## rangen = random.Random()
## nxcell = KeyTraits.maxKey1d/4
## assert nx < nxcell
## ncell = nxcell**2
## dxcell = 1.0/nxcell
## dycell = 1.0/nxcell
## occupiedCells = set()
## for i in xrange(n):
##     i = rangen.randint(0, ncell)
##     while i in occupiedCells:
##         i = rangen.randint(0, ncell)
##     ix = i % nxcell
##     iy = i // nxcell
##     generators.append(Vector((ix + 0.5)*dxcell, (iy + 0.5)*dycell))
##     occupiedCells.add(i)
## assert len(occupiedCells) == n
## nodes = createNodes(generators)
## print "Generating random mesh."
## mesh = PolygonalMesh(generators, Vector(0,0), Vector(1,1))
## print "Writing..."
## siloMeshDump("testPolygonalRandom", mesh, nodeLists = [nodes])
## print "Done."
## del nodes

## generators = vector_of_Vector()
## print "Creating generators for cylindrical mesh."
## dr = 1.0/nx
## for ir in xrange(nx):
##     ri = (ir + 0.5)*dr
##     ntheta = max(1, int(0.5*pi*ri / dr + 0.5))
##     dtheta = 0.5*pi/ntheta
##     for itheta in xrange(ntheta):
##         theta = (itheta + 0.5)*dtheta
##         generators.append(Vector(ri*cos(theta), ri*sin(theta)))
## print "Making nodes."
## nodes, db = createNodes(generators, dr)
## print "Generating cylindrical mesh."

# firstMom = VectorField("first moment", nodes)
# cm = db.connectivityMap()
# pos = nodes.positions()
# H = nodes.Hfield()
# for i in xrange(nodes.numInternalNodes):
#     wsum = 0.0
#     allneighbors = cm.connectivityForNode(nodes, i)
#     neighbors = allneighbors[0]
#     for j in neighbors:
#         etai = H[i]*(pos[j] - pos[i])
#         wi = W(etai, H[i])
#         wsum += wi
#         firstMom[i] += wi * etai
#     firstMom[i] /= wsum

# nPerh = nodes.nodesPerSmoothingScale
# void = makeVoidNodeList("void", 
#                         hmin = nodes.hmin,
#                         hmax = nodes.hmax,
#                         hminratio = nodes.hminratio,
#                         nPerh = nPerh)

# vpos = void.positions()
# vH = void.Hfield()
# for i in xrange(nodes.numInternalNodes):
#     fmi = firstMom[i]
#     if fmi.magnitude()/2.01 > 0.2:
#         fmhat = fmi.unitVector()
#         void.numInternalNodes += 1
#         j = void.numInternalNodes - 1
#         vpos[j] = pos[i] - H[i].Inverse()*fmhat/nPerh
#         vH[j] = H[i]
# print "Created %i void nodes." % void.numInternalNodes

## mesh, void = generatePolygonalMesh([nodes],
##                                    xmin = Vector(0.0, 0.0),
##                                    xmax = Vector(5.0, 5.0),
##                                    generateVoid = False,
##                                    generateParallelConnectivity = False,
##                                    removeBoundaryZones = True)


## # Compute the moments.
## print "Computing moments."
## nodeLists = vector_of_NodeList()
## for ns in [nodes, void]:
##     nodeLists.append(ns)
##     ns.neighbor().updateNodes()
## mom0 = ScalarFieldList(FieldListBase.Copy)
## mom1 = VectorFieldList(FieldListBase.Copy)
## zerothAndFirstNodalMoments(nodeLists, W, True, mom0, mom1)

## print "Writing..."
## siloMeshDump("testPolygonalCylindrical", mesh,
##              nodeLists = [nodes, void],
##              scalarFields = list(mom0),
##              vectorFields = list(mom1))
## print "Done."
