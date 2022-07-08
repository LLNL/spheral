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
  CXXTests
  Damage
  DataBase
  DataOutput
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

#spheral_add_cxx_library(CXX)
get_property(SPHERAL_CXX_LIBS GLOBAL PROPERTY SPHERAL_CXX_LIBS)
message("${SPHERAL_CXX_LIBS}")
message("CHECK")
blt_add_library(NAME        Spheral_CXX
                SOURCES     spheralCXX.cc
                DEPENDS_ON  ${SPHERAL_CXX_LIBS} ${SPHERAL_CXX_DEPENDS}
                SHARED      TRUE
                )

if(ENABLE_CUDA)
  set_target_properties(Spheral_CXX PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()

# Install Spheral C++ target and set it as an exportable CMake target
install(TARGETS             Spheral_CXX
        DESTINATION         lib
        EXPORT              ${PROJECT_NAME}-targets
        )

# Set the r-path of the C++ lib such that it is independent of the build dir when installed
set_target_properties(Spheral_CXX PROPERTIES
                    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${qhull_DIR}/lib;${conduit_DIR}/lib;${axom_DIR}/lib;${boost_DIR}/lib"
                    )

include_directories(${extra_packages_DIR})
foreach(e_package ${extra_packages})
  add_subdirectory(${extra_packages_DIR}/${e_package} ${CMAKE_CURRENT_BINARY_DIR}/${e_package})
endforeach()

install(EXPORT ${PROJECT_NAME}-targets DESTINATION Spheral/lib/cmake)
