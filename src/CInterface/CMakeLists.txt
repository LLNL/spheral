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

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_CInterface)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
