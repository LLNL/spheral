# CMakeLists to build the Spheral library.

cmake_minimum_required(VERSION 3.10)
project(spheral LANGUAGES CXX Fortran)

include(ExternalProject)

################################
# Configure CMake
################################
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -w")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_EXPORT_COMPILE_COMMANDS On)

################################
# Configure and Include blt
################################
set(ENABLE_MPI ON CACHE BOOL "")
set(ENABLE_OPENMP ON CACHE BOOL "")

if(NOT SPHERAL_BLT_DIR) 
  set (SPHERAL_BLT_REL_DIR "${PROJECT_SOURCE_DIR}/cmake/blt" CACHE PATH "")
  get_filename_component(SPHERAL_BLT_DIR "${SPHERAL_BLT_REL_DIR}" ABSOLUTE)
endif()

if (NOT EXISTS "${SPHERAL_BLT_DIR}/SetupBLT.cmake")
    message(FATAL_ERROR 
            "${SPHERAL_BLT_DIR} is not present.\n"
            "call cmake with -DSPHERAL_BLT_DIR=/your/installation/of/blt\n")
endif()

include(${SPHERAL_BLT_DIR}/SetupBLT.cmake)

################################
# Include standard build system logic and options / definitions
################################
set(ENABLE_CXXONLY OFF CACHE BOOL "enable C++ only build without python bindings")
set(ENABLE_2D ON CACHE BOOL "enable 2d")
set(ENABLE_3D ON CACHE BOOL "enable 3d")
set(ENABLE_INSTANTIATIONS ON CACHE BOOL "enable instantiations")
set(ENABLE_TIMER OFF CACHE BOOL "enable timer")

option(ENABLE_STATIC_CXXONLY "build only static libs" OFF)
if(ENABLE_STATIC_CXXONLY)
  set(ENABLE_CXXONLY ON)
endif()

if(ENABLE_MPI)
  set(BLT_MPI_COMPILE_FLAGS -DUSE_MPI -DMPICH_SKIP_MPICXX -ULAM_WANT_MPI2CPP -DOMPI_SKIP_MPICXX)
  list(APPEND spheral_blt_depends mpi)
endif()

if(ENABLE_OPENMP)
  list(APPEND spheral_blt_depends openmp)
endif()


################################
# Install / Locate third party libraries
################################
set(SPHERAL_INSTALL_DIR "" CACHE STRING "Directory to install Spheral TPLs and/or Spheral libs.")
include(cmake/InstallTPLs.cmake)

if(ENABLE_CXXONLY)
  set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/Spheral)
  if(SPHERAL_INSTALL_DIR)
    set(CMAKE_INSTALL_PREFIX ${SPHERAL_INSTALL_DIR})
  endif()
else()
  set(CMAKE_INSTALL_PREFIX ${PYTHON_SITE_PACKAGE_DIR})
endif()


include(cmake/CMakeDefinitions.cmake)


################################
# Set full rpath information by default
################################
# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_subdirectory(src)
