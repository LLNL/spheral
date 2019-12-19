###############################################################################
#
# Setup Eigen
# This file defines:
#  EIGEN_FOUND - If Eigen was found
#  EIGEN_INCLUDE_DIRS - The Eigen include directories
#  EIGEN_LIBRARY - The Eigen library

# first Check for EIGEN_DIR

if(NOT EIGEN_DIR)
    MESSAGE(FATAL_ERROR "Could not find EIGEN. EIGEN support needs explicit EIGEN_DIR")
endif()

#find includes
find_path( EIGEN_INCLUDE_DIRS Eigen
           PATHS  ${EIGEN_DIR}/include
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set EIGEN_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(EIGEN  DEFAULT_MSG
                                  EIGEN_INCLUDE_DIRS )
