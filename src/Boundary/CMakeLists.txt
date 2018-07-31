set(Boundary_sources2d
    AxisBoundaryRZ.cc
    BoundaryInst.cc.2.cc
    CRKSPHVoidBoundaryInst.cc.2.cc
    ConstantBoundaryInst.cc.2.cc
    ConstantRVelocityBoundaryInst.cc.2.cc
    ConstantVelocityBoundaryInst.cc.2.cc
    ConstantXVelocityBoundaryInst.cc.2.cc
    ConstantYVelocityBoundaryInst.cc.2.cc
    PeriodicBoundaryInst.cc.2.cc
    PlanarBoundaryInst.cc.2.cc
    ReflectingBoundaryInst.cc.2.cc
    RigidBoundaryInst.cc.2.cc
    )

set(Boundary_sources3d
    AxialSymmetryBoundary.cc
    BoundaryInst.cc.3.cc
    CRKSPHVoidBoundaryInst.cc.3.cc
    ConstantBoundaryInst.cc.3.cc
    ConstantRVelocityBoundaryInst.cc.3.cc
    ConstantVelocityBoundaryInst.cc.3.cc
    ConstantXVelocityBoundaryInst.cc.3.cc
    ConstantYVelocityBoundaryInst.cc.3.cc
    ConstantZVelocityBoundaryInst.cc.3.cc
    CylindricalBoundary.cc
    PeriodicBoundaryInst.cc.3.cc
    PlanarBoundaryInst.cc.3.cc
    ReflectingBoundaryInst.cc.3.cc
    RigidBoundaryInst.cc.3.cc
    SphericalBoundary.cc
    )

set(Boundary_headers
    AxialSymmetryBoundary.hh
    AxisBoundaryRZ.hh
    Boundary.hh
    BoundaryInline.hh
    CRKSPHVoidBoundary.hh
    ConstantBoundary.hh
    ConstantBoundaryInline.hh
    ConstantRVelocityBoundary.hh
    ConstantVelocityBoundary.hh
    ConstantVelocityBoundaryInline.hh
    ConstantXVelocityBoundary.hh
    ConstantYVelocityBoundary.hh
    ConstantZVelocityBoundary.hh
    CylindricalBoundary.hh
    CylindricalBoundaryInline.hh
    PeriodicBoundary.hh
    PeriodicPlanarBoundary.hh
    PlanarBoundary.hh
    PlanarBoundaryInline.hh
    PlanarSymmetryBoundary.hh
    ReflectingBoundary.hh
    ReflectingBoundaryInline.hh
    RigidBoundary.hh
    RigidBoundaryInline.hh
    SphericalBoundary.hh
    SphericalBoundaryInline.hh
    ThreeOneDimensionBoundary.hh
    mapPositionThroughPlanes.hh
    )

INCLUDE_DIRECTORIES(../)

set(spheral_depends BOOST Eigen qhull silo polytope mpi)

blt_add_library(NAME        Spheral_Boundary
                SOURCES     ${Boundary_sources2d} ${Boundary_sources3d}
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
