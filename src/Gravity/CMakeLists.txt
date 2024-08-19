include_directories(.)

set(Gravity_sources
    ApproximatePolyhedralGravityModel.cc
    )

# Add explicit instantiations
list(APPEND Gravity_sources
     TreeInst.cc
     TreeGravityInst.cc
     PolyGravityInst.cc
     NBodyGravityInst.cc
     )

set(Gravity_headers
  CompatibleGravitationalVelocityPolicy.hh
  NBodyGravity.hh
  Tree.hh
  TreeGravity.hh
  TreeGravityInline.hh
  ApproximatePolyhedralGravityModel.hh
  ApproximatePolyhedralGravityModelInline.hh
  PolyGravity.hh
  )

spheral_add_obj_library(Gravity SPHERAL_OBJ_LIBS)

