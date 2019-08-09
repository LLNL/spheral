set(Boundary_inst
    Boundary
    CRKSPHVoidBoundary
    PeriodicBoundary
    PlanarBoundary
    ReflectingBoundary
    RigidBoundary
   )


set(Boundary_sources )

if(ENABLE_2D)
  list(APPEND Boundary_sources
       AxisBoundaryRZ.cc
      )
endif()

if(ENABLE_3D)
  list(APPEND Boundary_sources
       AxialSymmetryBoundary.cc
      )
endif()

instantiate(Boundary_inst Boundary_sources)

set(Boundary_headers
    AxialSymmetryBoundary.hh
    AxisBoundaryRZ.hh
    Boundary.hh
    BoundaryInline.hh
    CRKSPHVoidBoundary.hh
    PeriodicBoundary.hh
    PeriodicPlanarBoundary.hh
    PlanarBoundary.hh
    PlanarBoundaryInline.hh
    PlanarSymmetryBoundary.hh
    ReflectingBoundary.hh
    ReflectingBoundaryInline.hh
    RigidBoundary.hh
    RigidBoundaryInline.hh
    mapPositionThroughPlanes.hh
    )

blt_add_library(NAME        Spheral_Boundary
                SOURCES     ${Boundary_sources}
                HEADERS     ${Boundary_headers}
                DEPENDS_ON  ${spheral_depends}
                )

install(TARGETS              Spheral_Boundary
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Boundary_headers}
        DESTINATION include/Boundary
        )
