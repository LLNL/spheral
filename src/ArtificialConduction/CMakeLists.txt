include_directories(.)

set(ArtificialConduction_sources
    ArtificialConduction.cc
    ArtificialConductionPolicy.cc
    )

set(ArtificialConduction_headers
    ArtificialConduction.hh
    ArtificialConductionPolicy.hh
    ) 

spheral_add_obj_library(ArtificialConduction SPHERAL_OBJ_LIBS)

