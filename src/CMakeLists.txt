cmake_minimum_required(VERSION 3.10)
project(Spheral)

#--------------------------------------------------------------------------------
#   Start of BLT configuration.
#--------------------------------------------------------------------------------

if (DEFINED BLT_SOURCE_DIR)
  if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
    message(FATAL_ERROR "Given BLT_SOURCE_DIR does not contain SetupBLT.cmake")
  endif()
else()
  # Use internal BLT if no BLT_SOURCE_DIR is given
  set(BLT_SOURCE_DIR "${PROJECT_SOURCE_DIR}/thirdPartyLibs/blt" CACHE PATH "")
  if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
    message(FATAL_ERROR
      "The BLT git submodule is not present in ${PROJECT_SOURCE_DIR}\n"
      "Either run the following two commands in your git repository: \n"
      "    git submodule init\n"
      "    git submodule update\n"
      "Or add -DBLT_SOURCE_DIR=/path/to/blt to your CMake command." )
  endif()
endif()

include(${BLT_SOURCE_DIR}/SetupBLT.cmake)

#--------------------------------------------------------------------------------
#   End of BLT configuration.
#--------------------------------------------------------------------------------

# Finding included packages.
find_package(GSL REQUIRED)
include(CMakeBasics.cmake)

# Add package subdirectories.
set(_packages
    ArtificialConduction
    Geometry
    Neighbor
    Utilities
    )

foreach(_package ${_packages})
  add_subdirectory(${_package})
endforeach()
