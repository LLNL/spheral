include_directories(.)

set(ArtificialConduction_sources )

# Add explicit instantiations
list(APPEND ArtificialConduction_sources
     ArtificialConductionInst.cc
     ArtificialConductionPolicyInst.cc
    )

set(ArtificialConduction_headers
    ArtificialConduction.hh
    ArtificialConductionPolicy.hh
    ) 

spheral_add_obj_library(ArtificialConduction SPHERAL_OBJ_LIBS)

