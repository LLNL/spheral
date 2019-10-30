set(Boundary_inst
  Boundary
  ConstantBoundary
  ConstantRVelocityBoundary
  ConstantVelocityBoundary
  ConstantXVelocityBoundary
  ConstantYVelocityBoundary
  CRKSPHVoidBoundary
  PeriodicBoundary
  PlanarBoundary
  ReflectingBoundary
  RigidBoundary
  )

set(Boundary_sources
  SphericalBoundary.cc
  #ThreeOneDimensionBoundary.cc
  ConstantZVelocityBoundaryInst.cc
  )
instantiate(Boundary_inst Boundary_sources)

set(Boundary_headers
  AxisBoundaryRZ.hh
  Boundary.hh
  BoundaryInline.hh
  ConstantBoundary.hh
  ConstantBoundaryInline.hh
  ConstantRVelocityBoundary.hh
  ConstantVelocityBoundary.hh
  ConstantVelocityBoundaryInline.hh
  ConstantXVelocityBoundary.hh
  ConstantYVelocityBoundary.hh
  ConstantZVelocityBoundary.hh
  CRKSPHVoidBoundary.hh
  CylindricalBoundary.hh
  CylindricalBoundaryInline.hh
  makefile.in
  mapPositionThroughPlanes.hh
  PeriodicBoundary.hh
  PeriodicPlanarBoundary.hh
  PlanarBoundary.hh
  PlanarBoundaryInline.hh
  PlanarSymmetryBoundary.hh
  ReflectingBoundary.hh
  ReflectingBoundaryInline.hh
  RigidBoundary.hh
  SphericalBoundary.hh
  SphericalBoundaryInline.hh
  #ThreeOneDimensionBoundary.hh
  )

set(Boundary_includes
  ${SPHERAL_SOURCE_DIR}
  ${SPHERAL_SOURCE_DIR}/Boundary
  ${EIGEN_LIB_LOCATION}
  ${PYBIND11_LIB_LOCATION}
  )


blt_add_library(NAME        Spheral_Boundary
                HEADERS     ${Boundary_headers}
                SOURCES     ${Boundary_sources}
                DEPENDS_ON  ${spheral_depends}
                INCLUDES    ${Boundary_includes}
                )

install(TARGETS             Spheral_Boundary
        EXPORT              ${PROJECT_NAME}-targets
        DESTINATION         lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Boundary_headers}
        DESTINATION include/Boundary
        )

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_Boundary)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
