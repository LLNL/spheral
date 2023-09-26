################################
# Build spheral packages 
################################
include(${SPHERAL_CMAKE_MODULE_PATH}/spheral/InstantiateCXX.cmake)
include(${SPHERAL_CMAKE_MODULE_PATH}/spheral/SpheralAddLibs.cmake)
include(${SPHERAL_CMAKE_MODULE_PATH}/spheral/SpheralInstallPythonFiles.cmake)

include_directories(.)
include_directories(${PROJECT_BINARY_DIR}/src)

set(SPHERAL_PYTHON_INSTALL ${PROJECT_BINARY_DIR}/lib)

# These are packages that must be made into targets and added as dependencies to SpheralCXX
list(APPEND _target_packages
  ArtificialViscosity
  ArtificialConduction
  Boundary
  CRKSPH
  Damage
  DataBase
  DataOutput
  DEM
  ExternalForce
  Field
  FieldOperations
  FileIO
  FSISPH
  Geometry
  Gravity
  Hydro
  Integrator
  Kernel
  KernelIntegrator
  Material
  Mesh
  Neighbor
  NodeGenerators
  NodeList
  Physics
  RK
  GSPH
  SolidMaterial
  SPH
  Strength
  SVPH
  Utilities
  )

# These are directories that should be added but not made into targets
list(APPEND _dir_packages
  CXXTests)

if(ENABLE_MPI)
  # Add Distributed target
  list(APPEND _target_packages Distributed)
endif()

if(NOT ENABLE_CXXONLY)
  list(APPEND _dir_packages SimulationControl PYB11)
endif()

configure_file(config.hh.in
  ${PROJECT_BINARY_DIR}/src/config.hh)

foreach(_package ${_dir_packages})
  add_subdirectory(${_package})
endforeach()

# Create a list of internal targets needed as dependencies for SpheralCXX
set(cxx_obj_list )
foreach(_package ${_target_packages})
  list(APPEND cxx_obj_list "Spheral_${_package}")
  add_subdirectory(${_package})
endforeach()
# Add target nested in CXXTests
list(APPEND cxx_obj_list "Spheral_DeviceTestLib")

set(CXX_sources spheralCXX.cc)
# Must use quotes when passing lists as inputs for functions
spheral_add_cxx_library(CXX "${cxx_obj_list}")

# This calls LLNLSpheralInstallObjs.cmake
if(EXISTS ${EXTERNAL_SPHERAL_OBJ_CMAKE})
  include(${EXTERNAL_SPHERAL_OBJ_CMAKE})
endif()
