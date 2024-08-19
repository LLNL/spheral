include_directories(.)

set(Mesh_sources
    findMatchingVertex.cc
    LineEdge.cc
    LineFace.cc
    LineMesh.cc
    LineZone.cc
    MeshConstructionUtilities.cc
    PolygonalEdge.cc
    PolygonalFace.cc
    PolygonalMesh.cc
    PolygonalZone.cc
    PolyhedralEdge.cc
    PolyhedralFace.cc
    PolyhedralMesh.cc
    PolyhedralZone.cc
    computeElementOwner.cc
    )

# Add explicit instantiations
list(APPEND Mesh_sources
     MeshPolicyInst.cc
     computeGeneratorsInst.cc
     findMatchingVertexInst.cc
     generateMeshInst.cc
     copy2polytopeInst.cc
     )

set(Mesh_headers
    CounterClockwiseCompareElements.hh
    Edge.hh
    EdgeInline.hh
    Face.hh
    FaceInline.hh
    Mesh.hh
    MeshConstructionUtilities.hh
    MeshInline.hh
    MeshPolicy.hh
    Node.hh
    NodeInline.hh
    Zone.hh
    ZoneInline.hh
    computeElementOwner.hh
    computeGenerators.hh
    copy2polytope.hh
    findMatchingVertex.hh
    generateMesh.hh
    )

spheral_install_python_files(
  generateMesh.py
  siloMeshDump.py
  )

spheral_add_obj_library(Mesh SPHERAL_OBJ_LIBS)

