include_directories(.)
set(Mesh_inst
    MeshPolicy
    computeGenerators
    findMatchingVertex
    generateMesh
    copy2polytope
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

set(Mesh_python_targets
  generateMesh.py
  siloMeshDump.py)

blt_add_library(NAME        Spheral_Mesh
                SOURCES     ${Mesh_sources}
                HEADERS     ${Mesh_headers}
                DEPENDS_ON  ${spheral_depends}
                SHARED      TRUE
                )

install(TARGETS              Spheral_Mesh
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Mesh_headers}
        DESTINATION include/Mesh
        )

if(NOT ENABLE_CXXONLY)
  install(FILES       ${Mesh_python_targets}
          DESTINATION .
          )
endif()

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_Mesh)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
