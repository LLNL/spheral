import time
from Spheral import *

#-------------------------------------------------------------------------------
# Generic method for creating a mesh.
#-------------------------------------------------------------------------------
def genericGenerateMesh(nodeLists,
                        boundaries,
                        xmin,
                        xmax,
                        meshGhostNodes,
                        generateVoid,
                        generateParallelConnectivity,
                        removeBoundaryZones,
                        voidThreshold,
                        DataBase,
                        vector_of_NodeList,
                        vector_of_Boundary,
                        Mesh,
                        makeVoidNodeList,
                        generateMesh):

    # If needed look up boundary ranges.
    if xmin is None or xmax is None:
        db = DataBase()
        for nodes in nodeLists:
            db.appendNodeList(nodes)
        exec("xmin0, xmax0 = Vector%id(), Vector%id()" % (db.nDim, db.nDim))
        db.boundingBox(xmin0, xmax0, ghost=False)
        delta = 0.1*(xmax0 - xmin0)
        xmin0 -= delta
        xmax0 += delta
        if xmin is None:
            xmin = xmin0
        if xmax is None:
            xmax = xmax0
        del db

    # Clip the range by any boundaries.
    for b in boundaries:
        b.clip(xmin, xmax)
    print "New range : %s %s" % (xmin, xmax)

    nodeListsVec = vector_of_NodeList()
    for x in nodeLists:
        nodeListsVec.append(x)
    mesh = Mesh()
    voidNodes = makeVoidNodeList("void", 
                                 hmin = nodeLists[0].hmin,
                                 hmax = nodeLists[0].hmax,
                                 hminratio = nodeLists[0].hminratio,
                                 nPerh = nodeLists[0].nodesPerSmoothingScale)
    nodeListsVec.append(voidNodes)

    boundVec = vector_of_Boundary()
    for x in boundaries:
        boundVec.append(x)

    t0 = time.time()
    generateMesh(nodeListsVec,
                 boundVec,
                 xmin, xmax, 
                 meshGhostNodes,
                 generateVoid,
                 generateParallelConnectivity,
                 removeBoundaryZones,
                 voidThreshold,
                 mesh,
                 voidNodes)
    print "Required %g seconds to generate mesh." % (time.time() - t0)
    del nodeListsVec
    return mesh, voidNodes

#-------------------------------------------------------------------------------
# Generate a 1-D tessellation of the nodes in a set of NodeLists.
#-------------------------------------------------------------------------------
def generateLineMesh(nodeLists,
                     boundaries = [],
                     xmin = None,
                     xmax = None,
                     meshGhostNodes = False,
                     generateVoid = True,
                     generateParallelConnectivity = False,
                     removeBoundaryZones = False,
                     voidThreshold = 2.0):
    return genericGenerateMesh(nodeLists,
                               boundaries,
                               xmin, 
                               xmax,
                               meshGhostNodes,
                               generateVoid,
                               generateParallelConnectivity,
                               removeBoundaryZones,
                               voidThreshold,
                               DataBase1d,
                               vector_of_NodeList1d,
                               vector_of_Boundary1d,
                               LineMesh,
                               makeVoidNodeList1d,
                               generateMesh1d)

#-------------------------------------------------------------------------------
# Generate a 2-D tessellation of the nodes in a set of NodeLists.
#-------------------------------------------------------------------------------
def generatePolygonalMesh(nodeLists,
                          boundaries = [],
                          xmin = None,
                          xmax = None,
                          meshGhostNodes = False,
                          generateVoid = True,
                          generateParallelConnectivity = False,
                          removeBoundaryZones = False,
                          voidThreshold = 2.0):
    return genericGenerateMesh(nodeLists,
                               boundaries,
                               xmin, 
                               xmax,
                               meshGhostNodes,
                               generateVoid,
                               generateParallelConnectivity,
                               removeBoundaryZones,
                               voidThreshold,
                               DataBase2d,
                               vector_of_NodeList2d,
                               vector_of_Boundary2d,
                               PolygonalMesh,
                               makeVoidNodeList2d,
                               generateMesh2d)

#-------------------------------------------------------------------------------
# Generate a 3-D tessellation of the nodes in a set of NodeLists.
#-------------------------------------------------------------------------------
def generatePolyhedralMesh(nodeLists,
                           boundaries = [],
                           xmin = None,
                           xmax = None,
                           meshGhostNodes = False,
                           generateVoid = True,
                           generateParallelConnectivity = False,
                           removeBoundaryZones = False,
                           voidThreshold = 2.0):
    return genericGenerateMesh(nodeLists,
                               boundaries,
                               xmin, 
                               xmax,
                               meshGhostNodes,
                               generateVoid,
                               generateParallelConnectivity,
                               removeBoundaryZones,
                               voidThreshold,
                               DataBase3d,
                               vector_of_NodeList3d,
                               vector_of_Boundary3d,
                               PolyhedralMesh,
                               makeVoidNodeList3d,
                               generateMesh3d)
