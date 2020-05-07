####################
# Other config defaults
####################

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_definitions("-DDEBUG=1")
else()
    add_definitions("-DDEBUG=0")
endif()

if (ENABLE_BOUNDCHECKING)
  add_definitions(-D_GLIBCXX_DEBUG=1)
endif()

if(ENABLE_CXXONLY)
  add_definitions(-DCXXONLY=1)
endif()
add_definitions(-DUSE_TETGEN=0)
add_definitions(-DUSE_TRIANGLE=0)
add_definitions(-DNOPOLYTOPE=1)
add_definitions(-DSPHERAL1D=1)
if(ENABLE_2D)
  add_definitions(-DSPHERAL2D=1)
endif()
if(ENABLE_3D)
  add_definitions(-DNOR3D=1)
  add_definitions(-DSPHERAL3D=1)
endif()
if(ENABLE_TIMER)
  add_definitions(-DTIMER=1)
endif()

if (ENABLE_MPI)
    add_definitions(-DUSE_MPI=1)
endif()
