set(Mesh_inst
    MeshPolicy
    computeGenerators
    findMatchingVertex
    generateMesh
   )


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

instantiate(Mesh_inst Mesh_sources)

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

blt_add_library(NAME        Spheral_Mesh
                SOURCES     ${Mesh_sources}
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
