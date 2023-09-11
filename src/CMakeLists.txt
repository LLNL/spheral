################################
# Build spheral packages 
################################
include(${SPHERAL_CMAKE_MODULE_PATH}/spheral/InstantiateCXX.cmake)
include(${SPHERAL_CMAKE_MODULE_PATH}/spheral/SpheralAddLibs.cmake)
include(${SPHERAL_CMAKE_MODULE_PATH}/spheral/SpheralInstallPythonFiles.cmake)

include_directories(.)
include_directories(${PROJECT_BINARY_DIR}/src)

set(SPHERAL_PYTHON_INSTALL ${PROJECT_BINARY_DIR}/lib)

list(APPEND _packages 
  ArtificialViscosity
  ArtificialConduction
  Boundary
  CRKSPH
  Damage
  DataBase
  DataOutput
  DEM
  Distributed
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

if(NOT ENABLE_CXXONLY)
  list(APPEND _packages 
    SimulationControl
    PYB11
    )
endif()

configure_file(config.hh.in
  ${PROJECT_BINARY_DIR}/src/config.hh)
# CXXTests should be added but not made into a target
add_subdirectory(CXXTests)
foreach(_package ${_packages})
  add_subdirectory(${_package})
endforeach()

include_directories(${extra_packages_DIR})
foreach(e_package ${extra_packages})
  add_subdirectory(${extra_packages_DIR}/${e_package} ${CMAKE_CURRENT_BINARY_DIR}/${e_package})
endforeach()

set(CXX_sources spheralCXX.cc)
set(cxx_obj_list )
foreach(_package ${_packages})
  list(APPEND cxx_obj_list "Spheral_${_package}")
endforeach()
# Must use quotes when passing lists as inputs for functions
spheral_add_cxx_library(CXX "${cxx_obj_list}")

install(EXPORT spheral_cxx-targets DESTINATION lib/cmake)
