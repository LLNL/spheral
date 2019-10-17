cmake_minimum_required(VERSION 3.10)
project(Spheral)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")
# Find Python
option(USE_PYTHON "Generate Python interface" ON)
if (USE_PYTHON)
  set(Python_ADDITIONAL_VERSIONS 2.7)
  if (PYTHON_EXE)
    set(PYTHON_EXECUTABLE ${PYTHON_EXE})
    include(FindPythonVersion)
    find_python_version()
    find_package(PythonLibs)
    string(REPLACE "/bin/python${PYTHON_VERSION}" "" PYTHON_ROOT ${PYTHON_EXE})
    set(PYTHON_INCLUDE_DIR "${PYTHON_ROOT}/include/site/python${PYTHON_VERSION}")
    set(PYTHON_INCLUDE_DIRS "${PYTHON_ROOT}/include/python${PYTHON_VERSION}")
    set(PYTHON_LIB_DIR "${PYTHON_ROOT}/lib")
    set(PYTHON_LIBRARY "${PYTHON_LIB_DIR}/libpython.a")
  else()
    find_package(PythonLibs REQUIRED)
    find_package(PythonInterp REQUIRED)
    set(PYTHON_VERSION "${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}")
  endif()

  if (PYTHON_EXECUTABLE)
    message(STATUS "Looking for Python...")
    include(FindPythonModule)
    find_python_module(PYB11Generator REQUIRED)
    find_python_module(pybind11 REQUIRED)
    set(PYB11_DIR "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps")
    set(PYB11_MODULE_NAME "polytope")
    set(PY_PYB11 1)
    include_directories("${PYTHON_INCLUDE_DIR}")
    include_directories("${PYTHON_INCLUDE_DIRS}")
    if (PY_PYB11)
      set(HAVE_PYTHON ON)
      message(STATUS "Found Python version ${PYTHON_VERSION}")
      message(STATUS "Using Python executable ${PYTHON_EXECUTABLE}")
    else()
      set(HAVE_PYTHON OFF)
    endif()
  endif()

  if (HAVE_PYTHON)
    include(FindPythonSitePackages)
    find_python_site_packages()
  endif()
else()
  set(HAVE_PYTHON OFF)
endif()

#--------------------------------------------------------------------------------
#   Start of BLT configuration.
#--------------------------------------------------------------------------------

if (DEFINED BLT_SOURCE_DIR)
  if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
    message(FATAL_ERROR "Given BLT_SOURCE_DIR does not contain SetupBLT.cmake")
  endif()
else()
  # Use internal BLT if no BLT_SOURCE_DIR is given
  set(BLT_SOURCE_DIR "${PROJECT_SOURCE_DIR}/src/thirdPartyLibs/blt" CACHE PATH "")
  if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
    message(FATAL_ERROR
      "The BLT git submodule is not present in ${PROJECT_SOURCE_DIR}\n"
      "Either run the following two commands in your git repository: \n"
      "    git submodule init\n"
      "    git submodule update\n"
      "Or add -DBLT_SOURCE_DIR=/path/to/blt to your CMake command." )
  endif()
endif()

include(${BLT_SOURCE_DIR}/SetupBLT.cmake)

#--------------------------------------------------------------------------------
#   End of BLT configuration.
#--------------------------------------------------------------------------------

#TODO:remove this
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -w")

# Finding included packages.
find_package(GSL REQUIRED)
set(EIGEN_LIB_LOCATION "${PROJECT_SOURCE_DIR}/src/thirdPartyLibs/eigen/")
set(QHULL_LIB_LOCATION "${PROJECT_SOURCE_DIR}/src/thirdPartyLibs/qhull/src/")
include(CMakeBasics.cmake)

set(SPHERAL_SOURCE_DIR "${PROJECT_SOURCE_DIR}/src")
add_subdirectory(src)
