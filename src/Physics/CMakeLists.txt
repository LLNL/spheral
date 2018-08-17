set(Physics_sources2d
    GenericBodyForceInst.cc.2.cc
    GenericHydroInst.cc.2.cc
    PhysicsInst.cc.2.cc
    )

set(Physics_sources3d
    GenericBodyForceInst.cc.3.cc
    GenericHydroInst.cc.3.cc
    PhysicsInst.cc.3.cc
    )

set(Physics_headers
    GenericBodyForce.hh
    GenericHydro.hh
    GenericHydroInline.hh
    Physics.hh
    PhysicsInline.hh
    )

INCLUDE_DIRECTORIES(../)

set(spheral_depends BOOST Eigen qhull silo polytope mpi)

blt_add_library(NAME        Spheral_Physics
                SOURCES     ${Physics_sources2d} ${Physics_sources3d}
                HEADERS     ${Physics_headers}
                DEPENDS_ON  ${spheral_depends}
                )

install(TARGETS              Spheral_Physics
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Physics_headers}
        DESTINATION include/Physics
        )
