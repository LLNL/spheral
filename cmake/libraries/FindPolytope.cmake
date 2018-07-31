###############################################################################
#
# Setup Polytope
# This file defines:
#  POLYTOPE_FOUND - If Polytope was found
#  POLYTOPE_INCLUDE_DIRS - The Polytope include directories
#  POLYTOPE_LIBRARY - The Polytope library

# first Check for POLYTOPE_DIR

if(NOT POLYTOPE_DIR)
    MESSAGE(FATAL_ERROR "Could not find POLYTOPE. POLYTOPE support needs explicit POLYTOPE_DIR")
endif()

#find includes
find_path( POLYTOPE_INCLUDE_DIRS polytope
           PATHS  ${POLYTOPE_DIR}/include/
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

find_library( POLYTOPE_LIBRARY NAMES polytope voro_2d voro_3d
              PATHS ${POLYTOPE_DIR}/lib
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)


include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set POLYTOPE_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(POLYTOPE  DEFAULT_MSG
                                  POLYTOPE_INCLUDE_DIRS
                                  POLYTOPE_LIBRARY )
