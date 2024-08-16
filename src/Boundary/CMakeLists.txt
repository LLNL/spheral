include_directories(.)

set(Boundary_sources 
  FacetedVolumeBoundary.cc
  BoundaryInst.cc
  CRKSPHVoidBoundaryInst.cc
  PlanarBoundaryInst.cc
  ReflectingBoundaryInst.cc
  RigidBoundaryInst.cc
  PeriodicBoundaryInst.cc
  ConstantVelocityBoundaryInst.cc
  ConstantXVelocityBoundaryInst.cc
  ConstantRVelocityBoundaryInst.cc
  ConstantBoundaryInst.cc
  findNodesTouchingThroughPlanesInst.cc
  InflowOutflowBoundaryInst.cc
  )

if(ENABLE_1D)
  list(APPEND Boundary_sources
       SphericalOriginBoundary.cc
       )
endif()

if(ENABLE_2D)
  list(APPEND Boundary_sources
       AxisBoundaryRZ.cc
       )
endif()

if(ENABLE_2D OR ENABLE_3D)
  list(APPEND Boundary_sources
       ConstantYVelocityBoundaryInst.cc
       )
endif()

if(ENABLE_3D)
  list(APPEND Boundary_sources
      CylindricalBoundary.cc  
      SphericalBoundary.cc
      ConstantZVelocityBoundaryInst.cc
      )
endif()

set(Boundary_headers
    AxisBoundaryRZ.hh
    SphericalOriginBoundary.hh
    Boundary.hh
    BoundaryInline.hh
    CRKSPHVoidBoundary.hh
    ConstantBoundary.hh
    ConstantBoundaryInline.hh
    ConstantBoundaryUtilities.hh
    ConstantRVelocityBoundary.hh
    ConstantVelocityBoundary.hh
    ConstantVelocityBoundaryInline.hh
    ConstantXVelocityBoundary.hh
    ConstantYVelocityBoundary.hh
    ConstantZVelocityBoundary.hh
    CylindricalBoundary.hh
    CylindricalBoundaryInline.hh
    FacetedVolumeBoundary.hh
    FacetedVolumeBoundaryInline.hh
    InflowOutflowBoundary.hh
    InflowOutflowBoundaryInline.hh
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
    ThreeOneDimensionBoundary.hh
    findNodesTouchingThroughPlanes.hh
    mapPositionThroughPlanes.hh
    )

spheral_install_python_files(ExpandingDomain.py)

spheral_add_obj_library(Boundary SPHERAL_OBJ_LIBS)

