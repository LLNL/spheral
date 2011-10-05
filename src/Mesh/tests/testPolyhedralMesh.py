#ATS:test(SELF,        label="PolyhedralMesh serial unit tests")
#ATS:test(SELF, np=8,  label="PolyhedralMesh serial (8 proc) tests")
#ATS:test(SELF, np=64,  label="PolyhedralMesh serial (64 proc) tests")

from math import *
import unittest

from Spheral3d import *
from generateMesh import *
from SpheralTestUtilities import fuzzyEqual

#===============================================================================
# Load mpi, and figure out how may domains to set up, and which domain we are.
#===============================================================================
import mpi
rank = mpi.rank
numDomains = mpi.procs
nxproc = int(numDomains**(1.0/3.0) + 0.1)
assert nxproc*nxproc*nxproc == numDomains
nxyproc = nxproc*nxproc

#===============================================================================
# Create a global random number generator.
#===============================================================================
import random
rangen = random.Random()

#===============================================================================
# Some boundary conditions.
#===============================================================================
x0, x1 = 0.0, 1.0
y0, y1 = 0.0, 1.0
z0, z1 = 0.0, 1.0

nx, ny, nz = 4*nxproc, 4*nxproc, 4*nxproc
nxy = nx*ny
n = nx*ny*nz
nperdomain = n / numDomains
nxcell = KeyTraits.maxKey1d/4
nycell = nxcell
nzcell = nxcell
nxycell = nxcell*nycell
assert nx < nxcell
ncell = nxcell*nycell*nzcell
dxcell = (x1 - x0)/nxcell
dycell = (y1 - y0)/nycell
dzcell = (z1 - z0)/nzcell

xmin = Vector(x0, y0, z0)
xmax = Vector(x1, y1, z1)

xbc0 = ReflectingBoundary(Plane(Vector(x0, y0, z0), Vector(1.0, 0.0, 0.0)))
ybc0 = ReflectingBoundary(Plane(Vector(x0, y0, z0), Vector(0.0, 1.0, 0.0)))
zbc0 = ReflectingBoundary(Plane(Vector(x0, y0, z0), Vector(0.0, 0.0, 1.0)))
bclist = [xbc0, ybc0, zbc0]
if numDomains > 1:
    bclist.append(BoundingVolumeDistributedBoundary.instance())
bclist = []

#===============================================================================
# Iterate the H's to something reasonable.
#===============================================================================
def iterateThoseHs(nodes):
    db = DataBase()
    db.appendNodeList(nodes)
    for bc in bclist:
        bc.setAllGhostNodes(db)
        nodes.neighbor().updateNodes()
    for bc in bclist:
        bc.finalizeGhostBoundary()
    nodes.neighbor().updateNodes()
    vecbound = vector_of_Boundary()
    for bc in bclist:
        vecbound.append(bc)
    WT = TableKernel(BSplineKernel(), 1000)
    smooth = SPHSmoothingScale()
    iterateIdealH(db, vecbound, WT, smooth,
                  tolerance = 1.0e-4)
    return

#===============================================================================
# A counter to help in creating unique NodeList names.
#===============================================================================
itest = 0

#===============================================================================
# Test class for tests to apply to all meshes.
#===============================================================================
class PolyhedralMeshGenericTests:

    #---------------------------------------------------------------------------
    # Test numbers of elements.
    #---------------------------------------------------------------------------
    def testPolyhedralMeshNums0(self):
        mesh, void = generatePolyhedralMesh([self.nodes],
                                            xmin = xmin,
                                            xmax = xmax,
                                            generateVoid = False,
                                            generateParallelConnectivity = False)
        assert mesh.numZones == self.nodes.numInternalNodes

    #---------------------------------------------------------------------------
    # Test element IDs.
    #---------------------------------------------------------------------------
    def testPolyhedralMeshElementIDs(self):
        mesh, void = generatePolyhedralMesh([self.nodes],
                                            xmin = xmin,
                                            xmax = xmax,
                                            generateVoid = False,
                                            generateParallelConnectivity = False)
        for i in xrange(self.nodes.numInternalNodes):
            node = mesh.node(i)
            assert node.ID == i
        for i in xrange(mesh.numEdges):
            edge = mesh.edge(i)
            assert edge.ID == i
        for i in xrange(mesh.numFaces):
            face = mesh.face(i)
            assert face.ID == i
        for i in xrange(mesh.numZones):
            zone = mesh.zone(i)
            assert zone.ID == i
        return

    #---------------------------------------------------------------------------
    # Test that the zones in the mesh correspond to the correct seed nodes.
    #---------------------------------------------------------------------------
    def testPolyhedralMeshZoneOrder(self):
        mesh, void = generatePolyhedralMesh([self.nodes],
                                            xmin = xmin,
                                            xmax = xmax,
                                            generateVoid = False,
                                            generateParallelConnectivity = False)
        assert mesh.numZones >= self.nodes.numInternalNodes
        for i in xrange(self.nodes.numInternalNodes):
            zonehull = mesh.zone(i).convexHull()
            self.failUnless(zonehull.contains(self.pos[i]),
                            "Zone hull does not contain generator:  %s\n     vertices: %s\n   distance: %g" %
                            (self.pos[i],
                             [str(x) for x in zonehull.vertices()],
                             zonehull.distance(self.pos[i])))
        return

    #---------------------------------------------------------------------------
    # Test the minimum scale.
    #---------------------------------------------------------------------------
    def testPolyhedralMeshMinimumScale(self):
        mesh, void = generatePolyhedralMesh([self.nodes],
                                            xmin = xmin,
                                            xmax = xmax)
        self.failUnless(mesh.minimumScale <= self.dxmin, 
                        "Scales don't match:  %g %g" % (mesh.minimumScale, self.dxmin))
        return

    #---------------------------------------------------------------------------
    # Test the parallel domain info.
    #---------------------------------------------------------------------------
    def testPolyhedralMeshParallel(self):
        mesh, void = generatePolyhedralMesh([self.nodes],
                                            xmin = xmin,
                                            xmax = xmax,
                                            generateParallelConnectivity = True)

        neighborDomains = [int(x) for x in mesh.neighborDomains]
        sharedNodes = []
        for ll in mesh.sharedNodes:
            sharedNodes.append([int(x) for x in ll])
        assert len(neighborDomains) == len(mesh.sharedNodes)

##         # Check the correct domains are talking to each other.
##         nxproc = int(sqrt(numDomains))
##         assert nxproc*nxproc == numDomains
##         ixproc = rank % nxproc
##         iyproc = rank / nxproc
##         neighborDomainsAnswer = []
##         for iy in xrange(max(0, iyproc - 1), min(nxproc, iyproc + 2)):
##             for ix in xrange(max(0, ixproc - 1), min(nxproc, ixproc + 2)):
##                 if not (ix == ixproc and iy == iyproc):
##                     neighborDomainsAnswer.append(ix + iy*nxproc)
##         ok = mpi.allreduce((neighborDomains == neighborDomainsAnswer), mpi.MIN)
##         self.failUnless(ok, "Strange neighbor domains for %i : %s ?= %s" % (rank, neighborDomains, neighborDomainsAnswer))

        # Check that the communicated mesh nodes are consistent.
        boxInv = xmax - xmin
        boxInv = Vector(1.0/boxInv.x, 1.0/boxInv.y, 1.0/boxInv.z)
        for sendProc in xrange(numDomains):
            numChecks = mpi.bcast(len(neighborDomains), root=sendProc)
            assert mpi.allreduce(numChecks, mpi.MIN) == mpi.allreduce(numChecks, mpi.MAX)
            for k in xrange(numChecks):
                if rank == sendProc:
                    ksafe = k
                else:
                    ksafe = 0
                recvProc = mpi.bcast(neighborDomains[ksafe], root=sendProc)
                recvHashes = mpi.bcast([hashPosition(mesh.node(i).position(), xmin, xmax, boxInv) for i in sharedNodes[ksafe]], root=sendProc)
                recvPos = mpi.bcast([str(mesh.node(i).position()) for i in sharedNodes[ksafe]], root=sendProc)
                ok = True
                msg = ""
                if rank == recvProc:
                    assert sendProc in neighborDomains
                    kk = neighborDomains.index(sendProc)
                    assert kk < len(sharedNodes)
                    ok = ([hashPosition(mesh.node(i).position(), xmin, xmax, boxInv) for i in sharedNodes[kk]] == recvHashes)
                    msg = ("Shared node indicies don't match %i %i\n   %s != %s\n    %s\n    %s" %
                           (rank, sendProc, 
                            str([hashPosition(mesh.node(i).position(), xmin, xmax, boxInv) for i in sharedNodes[kk]]),
                            recvHashes,
                            [str(mesh.node(i).position()) for i in sharedNodes[kk]],
                            recvPos))
                self.failUnless(mpi.allreduce(ok, mpi.MIN), msg)

        return

    #---------------------------------------------------------------------------
    # Test each mesh point hashes uniquely.
    #---------------------------------------------------------------------------
    def testPolyhedralMeshHash(self):
        mesh, void = generatePolyhedralMesh([self.nodes],
                                            xmin = xmin,
                                            xmax = xmax)
        pos = [mesh.zone(i).position() for i in xrange(mesh.numZones)] + [mesh.node(i).position() for i in xrange(mesh.numNodes)]
        boxInv = xmax - xmin
        boxInv = Vector(1.0/boxInv.x, 1.0/boxInv.y, 1.0/boxInv.z)
        hashes = [hashPosition(x, xmin, xmax, boxInv) for x in pos]
        blarg = zip(hashes, pos)
        blarg.sort()
        for i in xrange(len(blarg) - 1):
            hash0 = blarg[i][0]
            hash1 = blarg[i+1][0]
            self.failIf(hash0 == hash1,
                        "%i: Non-unique hash:  %i %i %s %s %s %s %g" % (mpi.rank, i, mesh.numZones, str(hash0), str(hash1), blarg[i][1], blarg[i+1][1], (blarg[i][1] - blarg[i+1][1]).magnitude()))
        return

    #---------------------------------------------------------------------------
    # Test the zones of the nodes.
    #---------------------------------------------------------------------------
    def testPolyhedralMeshNodeZones(self):
        mesh, void = generatePolyhedralMesh([self.nodes],
                                            xmin = xmin,
                                            xmax = xmax)
        answer = {}
        for inode in xrange(mesh.numNodes):
            answer[inode] = set()
        for izone in xrange(mesh.numZones):
            nodeIDs = mesh.zone(izone).nodeIDs
            for inode in nodeIDs:
                answer[inode].add(izone)
        for inode in xrange(mesh.numNodes):
            zoneIDs = mesh.node(inode).zoneIDs
            for izone in zoneIDs:
                self.failUnless(izone in answer[inode] or izone == PolyhedralMesh.UNSETID,
                                "Missing zone %i for set in node %i: %s %s" %
                                (izone, inode, [x for x in zoneIDs], answer[inode]))

    #---------------------------------------------------------------------------
    # Test the opposite zones across faces.
    #---------------------------------------------------------------------------
    def testPolyhedralMeshOppZones(self):
        mesh, void = generatePolyhedralMesh([self.nodes],
                                            xmin = xmin,
                                            xmax = xmax)
        answer = [[] for i in xrange(mesh.numFaces)]
        for izone in xrange(mesh.numZones):
            faces = mesh.zone(izone).faceIDs
            for iface in faces:
                answer[iface].append(izone)

        for iface in xrange(mesh.numFaces):
            face = mesh.face(iface)
            zoneIDs = answer[iface]
            assert len(zoneIDs) in (1, 2)
            if len(zoneIDs) == 2:
                self.failUnless(face.oppositeZoneID(zoneIDs[0]) == zoneIDs[1],
                                "Bad opposites:  (%i, %i) != (%i, %i)" %
                                (zoneIDs[0], zoneIDs[1],
                                 face.oppositeZoneID(zoneIDs[0]), face.oppositeZoneID(zoneIDs[1])))
                self.failUnless(face.oppositeZoneID(zoneIDs[1]) == zoneIDs[0],
                                "Bad opposites:  (%i, %i) != (%i, %i)" %
                                (zoneIDs[0], zoneIDs[1],
                                 face.oppositeZoneID(zoneIDs[0]), face.oppositeZoneID(zoneIDs[1])))
            else:
                assert face.oppositeZoneID(zoneIDs[0]) == PolyhedralMesh.UNSETID

#===============================================================================
# Create a uniformly spaced nodes/mesh.
#===============================================================================
class UniformPolyhedralMeshTests(unittest.TestCase, PolyhedralMeshGenericTests):

    #---------------------------------------------------------------------------
    # Create the NodeList we'll use for generating the mesh.
    #---------------------------------------------------------------------------
    def setUp(self):
        global itest
        eos = GammaLawGasMKS(5.0/3.0, 1.0)
        self.nodes = makeFluidNodeList("test nodes %i" % itest, eos,
                                       numInternal = nperdomain,
                                       nPerh = 2.01,
                                       hmin = 1e-5,
                                       hmax = 0.3)
        itest += 1
        self.pos = self.nodes.positions()
        self.H = self.nodes.Hfield()

        # Generate positions and split them up between domains appropriately.
        dxproc = (x1 - x0)/nxproc
        dyproc = (y1 - y0)/nxproc
        dzproc = (z1 - z0)/nxproc

        ixproc = rank % nxproc
        iyproc = (rank % nxyproc) / nxproc
        izproc = rank / nxyproc
        xminproc = Vector(x0 + ixproc*dxproc, y0 + iyproc*dyproc, z0 + izproc*dzproc)
        xmaxproc = Vector(x0 + (ixproc + 1)*dxproc, y0 + (iyproc + 1)*dyproc, z0 + (izproc + 1)*dzproc)

        dxavg = (x1 - x0)/nx
        dyavg = (y1 - y0)/ny
        dzavg = (z1 - z0)/nz
        self.dxmin = dxavg
        xyznodes = [v for v in (Vector(x0 + (i % nx +         0.5)*dxavg,
                                       y0 + ((i % nxy) / nx + 0.5)*dyavg,
                                       z0 + (i / nxy +        0.5)*dzavg) for i in xrange(n))
                    if testPointInBox(v, xminproc, xmaxproc)]
        assert len(xyznodes) == nperdomain
        assert mpi.allreduce(len(xyznodes), mpi.SUM) == n
        sys.stderr.write("%i nodes in range [%s:%s]\n" % (len(xyznodes), xminproc, xmaxproc))

        # We now have the positions for each domain appropriately divided, so shuffle
        # the local positions.
        #random.shuffle(xyznodes)

        # Now we can set the node conditions.
        for i in xrange(nperdomain):
            self.pos[i] = xyznodes[i]
            self.H[i] = SymTensor(1.0/(2.0*dxavg), 0.0, 0.0,
                                  0.0, 1.0/(2.0*dyavg), 0.0,
                                  0.0, 0.0, 1.0/(2.0*dzavg))
        self.nodes.neighbor().updateNodes()

        # Fix up the H's.
        #iterateThoseHs(self.nodes)
        return

    #---------------------------------------------------------------------------
    # Standard destructor.
    #---------------------------------------------------------------------------
    def tearDown(self):
        del self.nodes
        return

#===============================================================================
# Create randomly spaced set of nodes in the unit square.
#===============================================================================
class RandomPolyhedralMeshTests(unittest.TestCase, PolyhedralMeshGenericTests):

    #---------------------------------------------------------------------------
    # Create the NodeList we'll use for generating the mesh.
    #---------------------------------------------------------------------------
    def setUp(self):
        global itest
        eos = GammaLawGasMKS(5.0/3.0, 1.0)
        self.nodes = makeFluidNodeList("test nodes %i" % itest, eos,
                                       numInternal = nperdomain,
                                       nPerh = 2.01,
                                       hmin = 1.0e-5,
                                       hmax = 0.3)
        itest += 1
        self.pos = self.nodes.positions()
        self.H = self.nodes.Hfield()

        # Figure out the domain bounding volumes.
        dxproc = (x1 - x0)/nxproc
        dyproc = (y1 - y0)/nxproc
        dzproc = (z1 - z0)/nxproc

        ixproc = rank % nxproc
        iyproc = (rank % nxyproc) / nxproc
        izproc = rank / nxyproc
        xminproc = Vector(x0 + ixproc*dxproc, y0 + iyproc*dyproc, z0 + izproc*dzproc)
        xmaxproc = Vector(x0 + (ixproc + 1)*dxproc, y0 + (iyproc + 1)*dyproc, z0 + (izproc + 1)*dzproc)

        # Randomly seed the generators.  We choose from random cells in order
        # to keep nodes from getting too close together.
        xyznodes_all = []
        occupiedCells = set()
        for k in xrange(n):
            i = rangen.randint(0, ncell)
            while i in occupiedCells:
                i = rangen.randint(0, ncell)
            ix = i % nxcell
            iy = (i % nxycell) / nxcell
            iz = i / nxycell
            xyznodes_all.append(Vector((ix + 0.5)*dxcell, (iy + 0.5)*dycell, (iz + 0.5)*dzcell))
            occupiedCells.add(i)
        assert len(occupiedCells) == n
        xyznodes_all = mpi.bcast(xyznodes_all)
        xyznodes = [v for v in xyznodes_all if testPointInBox(v, xminproc, xmaxproc)]
        dxavg = (x1 - x0)/nx
        dyavg = (y1 - y0)/ny
        dzavg = (z1 - z0)/nz
        self.dxmin = dxavg
        assert mpi.allreduce(len(xyznodes), mpi.SUM) == n

        # Now we can set the node conditions.
        self.nodes.numInternalNodes = len(xyznodes)
        for i in xrange(len(xyznodes)):
            self.pos[i] = xyznodes[i]
            self.H[i] = SymTensor(1.0/(2.0*dxavg), 0.0, 0.0,
                                  0.0, 1.0/(2.0*dyavg), 0.0,
                                  0.0, 0.0, 1.0/(2.0*dzavg))
        self.nodes.neighbor().updateNodes()

        # Fix up the H's.
        #iterateThoseHs(self.nodes)
        return

    #---------------------------------------------------------------------------
    # Standard destructor.
    #---------------------------------------------------------------------------
    def tearDown(self):
        del self.nodes
        return

#===============================================================================
# Run the tests
#===============================================================================
if __name__ == "__main__":
    unittest.main()
