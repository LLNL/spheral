include_directories(.)

set(Physics_sources )

# Add explicit instantiations
list(APPEND Physics_sources
     GenericBodyForceInst.cc
     GenericHydroInst.cc
     PhysicsInst.cc
     )

set(Physics_headers
    GenericBodyForce.hh
    GenericHydro.hh
    GenericHydroInline.hh
    Physics.hh
    PhysicsInline.hh
    )

spheral_add_obj_library(Physics SPHERAL_OBJ_LIBS)

