include_directories(.)
set(Gravity_inst
  CompatibleGravitationalVelocityPolicy
  NBodyGravity
  )

set(Gravity_sources
  TreeGravityInst.cc
  ApproximatePolyhedralGravityModel.cc
  )
instantiate(Gravity_inst Gravity_sources)

set(Gravity_headers
  CompatibleGravitationalVelocityPolicy.hh
  NBodyGravity.hh
  TreeGravity.hh
  TreeGravityInline.hh
  ApproximatePolyhedralGravityModel.hh
  ApproximatePolyhedralGravityModelInline.hh
  )


spheral_add_obj_library(Gravity  SPHERAL_OBJ_LIBS)

