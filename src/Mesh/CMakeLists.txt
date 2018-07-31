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

set(Mesh_sources2d
    MeshPolicyInst.cc.2.cc
    computeGeneratorsInst.cc.2.cc
    findMatchingVertexInst.cc.2.cc
    generateMeshInst.cc.2.cc
    )

set(Mesh_sources3d
    MeshPolicyInst.cc.3.cc
    computeGeneratorsInst.cc.3.cc
    findMatchingVertexInst.cc.3.cc
    generateMeshInst.cc.3.cc
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
    findMatchingVertex.hh
    generateMesh.hh
    )

INCLUDE_DIRECTORIES(../)

set(spheral_depends BOOST Eigen qhull silo polytope mpi)

blt_add_library(NAME        Spheral_Mesh
                SOURCES     ${Mesh_sources} ${Mesh_sources2d} ${Mesh_sources3d}
                HEADERS     ${Mesh_headers}
                DEPENDS_ON  ${spheral_depends}
                )

install(TARGETS              Spheral_Mesh
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Mesh_headers}
        DESTINATION include/Mesh
        )
