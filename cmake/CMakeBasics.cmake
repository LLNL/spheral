###############################################################################
# Include project's CMake files
###############################################################################
# Basic

# Libraries
include(cmake/SetupLibraries.cmake)

####################
# Other config defaults
####################

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_definitions("-DDEBUG=1")
else()
    add_definitions("-DDEBUG=0")
endif()

add_definitions(-DCXXONLY=1)
add_definitions(-DNOR3D=1)
add_definitions(-DUSE_TETGEN=0)
add_definitions(-DUSE_TRIANGLE=0)
add_definitions(-DNOPOLYTOPE=1)
add_definitions(-DSPHERAL2D=1)
add_definitions(-DSPHERAL3D=1)
add_definitions(-DTIMER=1)

add_compile_options(${SPHERAL_BGQ_FLAGS})

if (ENABLE_MPI)
    add_definitions(-DUSE_MPI=1)
endif()
