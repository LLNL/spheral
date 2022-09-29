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
  PolyGravity.hh
  )

# The PolyGravity object requires 3D, and optionally works in 2D
if(ENABLE_3D)
  list(APPEND Gravity_sources
       PolyGravity.cc)
endif()

spheral_add_obj_library(Gravity)

