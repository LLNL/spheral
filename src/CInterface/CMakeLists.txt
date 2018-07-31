set(CInterface_sources
    SpheralC.cc
    SpheralPseudoScript.cc
    )

set(CInterface_headers
    SpheralC.h
    SpheralPseudoScript.hh
    )

INCLUDE_DIRECTORIES(../)

set(spheral_depends BOOST Eigen qhull silo polytope mpi)

blt_add_library(NAME        Spheral_CInterface
                SOURCES     ${CInterface_sources}
                HEADERS     ${CInterface_headers}
                DEPENDS_ON  ${spheral_depends}
                )

install(TARGETS              Spheral_CInterface
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${CInterface_headers}
        DESTINATION include/CInterface
        )
