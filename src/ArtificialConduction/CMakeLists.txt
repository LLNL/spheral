include_directories(.)

set(ArtificialConduction_headers
    ArtificialConduction.hh
    ArtificialConductionDecl.hh
    ArtificialConductionImpl.hh
    ArtificialConductionPolicy.hh
    ArtificialConductionPolicyDecl.hh
    ArtificialConductionPolicyImpl.hh
    ) 

set(ArtificialConduction_sources )

if (ENABLE_INSTANTIATIONS)
  list(APPEND ArtificialConduction_sources
       ArtificialConductionInst.cc
       ArtificialConductionPolicyInst.cc
       )
endif ()

spheral_add_obj_library(ArtificialConduction SPHERAL_OBJ_LIBS)

