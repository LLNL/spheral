include_directories(.)
set(CInterface_sources
    SpheralC.cc
    SpheralPseudoScript.cc
    )

set(CInterface_headers
    SpheralC.h
    SpheralPseudoScript.hh
    )

spheral_add_cxx_library(CInterface)

install(TARGETS              Spheral_CInterface
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${CInterface_headers}
        DESTINATION include/CInterface
        )

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_CInterface)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
