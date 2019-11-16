# CMakeLists to build the Spheral library.
cmake_minimum_required(VERSION 3.4)
project(spheral LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -w")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")

include(cmake/InstallLibraries.cmake)

################################
# Include blt
################################
set (SPHERAL_BLT_REL_DIR "${PROJECT_SOURCE_DIR}/blt" CACHE PATH "")
get_filename_component(SPHERAL_BLT_DIR "${SPHERAL_BLT_REL_DIR}" ABSOLUTE)

if (NOT EXISTS "${PROJECT_SOURCE_DIR}/blt/SetupBLT.cmake")
    message(FATAL_ERROR 
            "${SPHERAL_BLT_DIR} is not present.\n"
            "call cmake with -DSPHERAL_BLT_DIR=/your/installation/of/blt\n")
endif()

include(${PROJECT_SOURCE_DIR}/blt/SetupBLT.cmake)

option(ENABLE_2D "enable 2d" ON)
option(ENABLE_3D "enable 3d" ON)
option(ENABLE_INSTANTIATIONS "enable instantiations" ON)
option(BUILD_PYTHON_TARGETS "enable build of spheral python tagets" ON)

################################
# Include standard build system logic and options
################################
include(cmake/CMakeBasics.cmake)

add_subdirectory(src)

install(EXPORT ${PROJECT_NAME}-targets DESTINATION lib/cmake)
