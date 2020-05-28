include_directories(.)
set(Gravity_inst
  CompatibleGravitationalVelocityPolicy
  NBodyGravity
  )

set(Gravity_sources
  #LaplacianMatrixFactory.cc
  TreeGravityInst.cc
  )
instantiate(Gravity_inst Gravity_sources)

set(Gravity_headers
  CompatibleGravitationalVelocityPolicy.hh
  #LaplacianMatrixFactory.hh
  NBodyGravity.hh
  TreeGravity.hh
  TreeGravityInline.hh
  )


spheral_add_cxx_library(Gravity)

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_Gravity)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
