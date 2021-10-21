# CMakeLists to build the Spheral library.

cmake_minimum_required(VERSION 3.10)
project(spheral LANGUAGES CXX Fortran)

set(SPHERAL_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

include(cmake/SetupSpheral.cmake)

if (ENABLE_TESTS)
  add_subdirectory(${SPHERAL_ROOT_DIR}/tests/Unit/CXXTests)

  install(DIRECTORY ${SPHERAL_ROOT_DIR}/tests/
        DESTINATION ${CMAKE_INSTALL_PREFIX}/tests
        PATTERN "runCXXTests.ats" EXCLUDE
       )
endif()

add_subdirectory(${SPHERAL_ROOT_DIR}/scripts)
