# CMakeLists to build the Spheral library.

cmake_minimum_required(VERSION 3.14)
include(cmake/SpheralVersion.cmake)
project(spheral LANGUAGES C CXX Fortran VERSION ${SPHERAL_VERSION})

set(SPHERAL_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE PATH "Path to Spheral source directory")
set(SPHERAL_TEST_INSTALL_PREFIX "tests")

include(cmake/SetupSpheral.cmake)

add_subdirectory(${SPHERAL_ROOT_DIR}/scripts)
