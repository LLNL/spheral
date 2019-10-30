set(Mesh_inst
  computeGenerators
  copy2polytope
  findMatchingVertex
  generateMesh
  MeshPolicy
  )

set(Mesh_sources
  computeElementOwner.cc
  findMatchingVertexTempl.cc
  LineEdge.cc
  LineFace.cc
  LineMesh.cc
  LineZone.cc
  Mesh.cc
  MeshConstructionUtilities.cc
  PolygonalEdge.cc
  PolygonalFace.cc
  PolygonalMesh.cc
  PolygonalZone.cc
  PolyhedralEdge.cc
  PolyhedralFace.cc
  PolyhedralMesh.cc
  PolyhedralZone.cc
  #siloMeshDump.cc
)
instantiate(Mesh_inst Mesh_sources)

set(Mesh_headers
  computeElementOwner.hh
  computeGenerators.hh
  copy2polytope.hh
  CounterClockwiseCompareElements.hh
  Edge.hh
  EdgeInline.hh
  Face.hh
  FaceInline.hh
  findMatchingVertex.hh
  generateMesh.hh
  Mesh.hh
  MeshConstructionUtilities.hh
  MeshInline.hh
  MeshPolicy.hh
  Node.hh
  NodeInline.hh
  Zone.hh
  ZoneInline.hh
)

set(Mesh_includes
  ${SPHERAL_SOURCE_DIR}
  ${SPHERAL_SOURCE_DIR}/Mesh
  ${EIGEN_LIB_LOCATION}
  ${PYBIND11_LIB_LOCATION}
  ${POLYTOPE_LIB_LOCATION}
  ${Boost_INCLUDE_DIRS}
  )


blt_add_library(NAME        Spheral_Mesh
                HEADERS     ${Mesh_headers}
                SOURCES     ${Mesh_sources}
                DEPENDS_ON  ${spheral_depends}
                INCLUDES    ${Mesh_includes}
                )

install(TARGETS             Spheral_Mesh
        EXPORT              ${PROJECT_NAME}-targets
        DESTINATION         lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Mesh_headers}
        DESTINATION include/Mesh
        )

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_Mesh)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
