# CMakeLists to build the Spheral library.
cmake_minimum_required(VERSION 3.4)
project(spheral LANGUAGES CXX)

################################
# Configure CMake
################################
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -w")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../cmake/Modules")


################################
# Configure and Include blt
################################
set(ENABLE_MPI ON CACHE BOOL "")

set (SPHERAL_BLT_REL_DIR "${PROJECT_SOURCE_DIR}/../blt" CACHE PATH "")
get_filename_component(SPHERAL_BLT_DIR "${SPHERAL_BLT_REL_DIR}" ABSOLUTE)

if (NOT EXISTS "${PROJECT_SOURCE_DIR}/../blt/SetupBLT.cmake")
    message(FATAL_ERROR 
            "${SPHERAL_BLT_DIR} is not present.\n"
            "call cmake with -DSPHERAL_BLT_DIR=/your/installation/of/blt\n")
endif()

include(${PROJECT_SOURCE_DIR}/../blt/SetupBLT.cmake)


################################
# Install / Locate third party libraries
################################
option(INSTALL_TPLS "download and install third party libs" ON)
include(../cmake/InstallLibraries.cmake)
set(spheral_depends BOOST Eigen qhull silo polytope pybind11)


################################
# Include standard build system logic and options / definitions
################################
option(ENABLE_CXXONLY "enable C++ only build without python bindings" OFF)
option(ENABLE_2D "enable 2d" ON)
option(ENABLE_3D "enable 3d" ON)
option(ENABLE_INSTANTIATIONS "enable instantiations" ON)

if(ENABLE_CXXONLY)
  set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/Spheral)
else()
  set(CMAKE_INSTALL_PREFIX ${PYTHON_SITE_PACKAGE_DIR}/Spheral)
endif()

if(ENABLE_MPI)
  set(BLT_MPI_COMPILE_FLAGS -DUSE_MPI -DMPICH_SKIP_MPICXX -ULAM_WANT_MPI2CPP -DOMPI_SKIP_MPICXX)
  list(APPEND spheral_depends mpi)
endif()

include(../cmake/CMakeBasics.cmake)


################################
# Build spheral packages 
################################
include_directories(.)
set(SPHERAL_PYTHON_DEPENDS )

set(SPHERAL_PYTHON_INSTALL ${PROJECT_BINARY_DIR}/lib)

set(_packages 
  ArtificialViscosity
  ArtificialConduction
  Boundary
  CRKSPH
  Damage
  DataBase
  DataOutput
  Distributed
  ExternalForce
  Field
  FieldOperations
  FileIO
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
  SolidMaterial
  SPH
  Strength
  SVPH
  Utilities
  )

if(NOT ENABLE_CXXONLY)
  list(APPEND _packages 
    PythonMPIInterfaces
    SimulationControl
    Pybind11Wraps
    )
endif()

if(NOT ENABLE_MPI)
  list(REMOVE_ITEM _packages
    Distributed
    )
endif()

foreach(_package ${_packages})
   add_subdirectory(${_package})
endforeach()

install(EXPORT ${PROJECT_NAME}-targets DESTINATION lib/cmake)
