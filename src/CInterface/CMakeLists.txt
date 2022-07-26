include_directories(.)
set(CInterface_sources
    SpheralC.cc
    SpheralPseudoScript.cc
    )

set(CInterface_headers
    SpheralC.h
    SpheralPseudoScript.hh
    )

spheral_add_obj_library(CInterface)
