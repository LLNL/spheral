include_directories(.)

set(ExternalForce_sources )

# Add explicit instantiations
list(APPEND ExternalForce_sources
     ConstantAccelerationInst.cc
     LinearAccelerationInst.cc
     NFWPotentialInst.cc
     PointPotentialInst.cc
     )

set(ExternalForce_headers
    ConstantAcceleration.hh
    ConstantAccelerationInline.hh
    LinearAcceleration.hh
    LinearAccelerationInline.hh
    NFWPotential.hh
    NFWPotentialInline.hh
    PointPotential.hh
    PointPotentialInline.hh
    )

set(ExternalForce_includes
   ${SPHERAL_SOURCE_DIR}
   ${SPHERAL_SOURCE_DIR}/ExternalForce
   ${EIGEN_LIB_LOCATION}
   )

spheral_add_obj_library(ExternalForce SPHERAL_OBJ_LIBS)
