###############################################################################
#
# Setup Pybind11
# This file defines:
#  PYBIND11_FOUND - If Pybind11 was found
#  PYBIND11_INCLUDE_DIRS - The Pybind11 include directories
#  PYBIND11_LIBRARY - The Pybind11 library

# first Check for PYBIND11_DIR

if(NOT PYBIND11_DIR)
    MESSAGE(FATAL_ERROR "Could not find PYBIND11. PYBIND11 support needs explicit PYBIND11_DIR")
endif()

#find includes
find_path( PYBIND11_INCLUDE_DIRS pybind11
           PATHS  ${PYBIND11_DIR}/include
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set PYBIND11_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(PYBIND11  DEFAULT_MSG
                                  PYBIND11_INCLUDE_DIRS )
