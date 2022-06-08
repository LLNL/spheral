################################
# Build spheral packages 
################################
include(${CMAKE_MODULE_PATH}/spheral/InstantiateCXX.cmake)
include(${CMAKE_MODULE_PATH}/spheral/SpheralAddLibs.cmake)
include(${CMAKE_MODULE_PATH}/spheral/SpheralInstallPythonFiles.cmake)

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
    Pybind11Wraps
    )
endif()

configure_file(config.hh.in
  ${PROJECT_BINARY_DIR}/src/config.hh)

foreach(_package ${_packages})
  add_subdirectory(${_package})
endforeach()

include_directories(${extra_packages_DIR})
foreach(e_package ${extra_packages})
  add_subdirectory(${extra_packages_DIR}/${e_package} ${CMAKE_CURRENT_BINARY_DIR}/${e_package})
endforeach()

install(EXPORT ${PROJECT_NAME}-targets DESTINATION Spheral/lib/cmake)
