###############################################################################
#
# Setup Qhull
# This file defines:
#  QHULL_FOUND - If Qhull was found
#  QHULL_INCLUDE_DIRS - The Qhull include directories
#  QHULL_LIBRARY - The Qhull library

# first Check for QHULL_DIR

if(NOT QHULL_DIR)
    MESSAGE(FATAL_ERROR "Could not find QHULL. QHULL support needs explicit QHULL_DIR")
endif()

#find includes
find_path( QHULL_INCLUDE_DIRS libqhull
           PATHS  ${QHULL_DIR}/include/
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

find_library( QHULL_LIBRARY NAMES qhullstatic
              PATHS ${QHULL_DIR}/lib
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)


include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set QHULL_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(QHULL  DEFAULT_MSG
                                  QHULL_INCLUDE_DIRS
                                  QHULL_LIBRARY )
