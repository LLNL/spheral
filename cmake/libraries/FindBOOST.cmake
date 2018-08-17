###############################################################################
#
# Setup Boost
# This file defines:
#  BOOST_FOUND - If Boost was found
#  BOOST_INCLUDE_DIRS - The Boost include directories
#  BOOST_LIBRARY - The Boost library

# first Check for BOOST_DIR

if(NOT BOOST_DIR)
    MESSAGE(FATAL_ERROR "Could not find BOOST. BOOST support needs explicit BOOST_DIR")
endif()

#find includes
find_path( BOOST_INCLUDE_DIRS boost
           PATHS  ${BOOST_DIR}/include
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set BOOST_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(BOOST  DEFAULT_MSG
                                  BOOST_INCLUDE_DIRS )
