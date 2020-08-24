include_directories(.)
set(Boundary_inst
    Boundary
    CRKSPHVoidBoundary
    PlanarBoundary
    ReflectingBoundary
    RigidBoundary
    PeriodicBoundary
    ConstantVelocityBoundary
    ConstantXVelocityBoundary
    ConstantYVelocityBoundary
    ConstantRVelocityBoundary
    ConstantBoundary
    findNodesTouchingThroughPlanes
    InflowOutflowBoundary
   )


set(Boundary_sources 
  FacetedVolumeBoundary.cc
       )

if(ENABLE_2D)
  list(APPEND Boundary_sources
       AxisBoundaryRZ.cc
      )
endif()

if(ENABLE_3D)
  list(APPEND Boundary_sources
      CylindricalBoundary.cc  
      SphericalBoundary.cc
      ConstantZVelocityBoundaryInst.cc
      )
endif()

instantiate(Boundary_inst Boundary_sources)

set(Boundary_headers
  #AxialSymmetryBoundary.hh
    AxisBoundaryRZ.hh
    Boundary.hh
    BoundaryInline.hh
    CRKSPHVoidBoundary.hh
    PlanarBoundary.hh
    PlanarBoundaryInline.hh
    PlanarSymmetryBoundary.hh
    ReflectingBoundary.hh
    ReflectingBoundaryInline.hh
    RigidBoundary.hh
    #RigidBoundaryInline.hh
    mapPositionThroughPlanes.hh
    )

spheral_add_cxx_library(Boundary)

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_Boundary)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
