message("-- C++ Compiler ID: ${CMAKE_CXX_COMPILER_ID}")

#-------------------------------------------------------------------------------
# Optionally suppress compiler warnings
#-------------------------------------------------------------------------------
option(ENABLE_WARNINGS "show compiler warnings" ON)
option(ENABLE_WARNINGS_AS_ERRORS "make warnings errors" OFF)

option(ENABLE_UNUSED_VARIABLE_WARNINGS "show unused variable compiler warnings" ON)
option(ENABLE_UNUSED_PARAMETER_WARNINGS "show unused parameter warnings" OFF)
option(ENABLE_MISSING_INCLUDE_DIR_WARNINGS "show unused parameter warnings" ON)


set(CXX_WARNING_FLAGS "")
if (ENABLE_WARNINGS)
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    list(APPEND CXX_WARNING_FLAGS -fdiagnostics-show-option -Wno-unused-command-line-argument -Wno-c++17-extensions)
    if(CMAKE_CXX_COMPILER_VERSION GREATER_EQUAL 18.0.0)
      list(APPEND CXX_WARNING_FLAGS -Wno-enum-constexpr-conversion)
    endif()
  endif()
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
endif()
message("-- Compiler warnings ${ENABLE_WARNINGS}")

if (ENABLE_WARNINGS_AS_ERRORS)
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    list(APPEND CXX_WARNING_FLAGS /W4 /WX)
  else()
    list(APPEND CXX_WARNING_FLAGS -Wall -Wextra -pedantic -Werror -Wl,--fatal-warnings)
  endif()
  message("-- Treating warnings as errors with compile flags ${CXX_WARNING_FLAGS}")
endif()


if (NOT ENABLE_UNUSED_VARIABLE_WARNINGS)
  list(APPEND CXX_WARNING_FLAGS -Wno-unused-variable)
endif()
message("-- Compiler unused variable warnings ${ENABLE_UNUSED_VARIABLE_WARNINGS}")


if (NOT ENABLE_UNUSED_PARAMETER_WARNINGS)
  list(APPEND CXX_WARNING_FLAGS -Wno-unused-parameter)
endif()
message("-- Compiler unused parameter warnings ${ENABLE_UNUSED_PARAMETER_WARNINGS}")


if (NOT ENABLE_MISSING_INCLUDE_DIR_WARNINGS)
  list(APPEND CXX_WARNING_FLAGS -Wno-missing-include-dirs)
endif()
message("-- Compiler missing include dir warnings ${ENABLE_MISSING_INCLUDE_DIR_WARNINGS}")

set(CUDA_WARNING_FLAGS -Xcudafe=\"--diag_suppress=esa_on_defaulted_function_ignored\")

add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:${CXX_WARNING_FLAGS}>")
add_compile_options("$<$<COMPILE_LANGUAGE:CUDA>:${CUDA_WARNING_FLAGS}>")
message("-- using warning flags ${CXX_WARNING_FLAGS}")

# We build some Fortran code from outside sources (like the Helmholtz EOS) that
# cause building errors if the compiler is too picky...
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-missing-include-dirs")
message("-- Fortran flags: ${CMAKE_Fortran_FLAGS}")

#-------------------------------------------------------------------------------
# PYB11 Target Flags
#-------------------------------------------------------------------------------
set(SPHERAL_PYB11_TARGET_FLAGS
  -O1
  -Wno-unused-local-typedefs 
  -Wno-overloaded-virtual)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  list(APPEND SPHERAL_PYB11_TARGET_FLAGS
    -Wno-self-assign-overloaded 
    -Wno-inconsistent-missing-override
    -Wno-delete-non-abstract-non-virtual-dtor
    -Wno-delete-abstract-non-virtual-dtor)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  list(APPEND SPHERAL_PYB11_TARGET_FLAGS
    -Wno-pedantic
    -fno-var-tracking-assignments)
endif()

#-------------------------------------------------------------------------------
# Compiler specific flags
#-------------------------------------------------------------------------------
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS -wd11074,11076,654)
  set(SPHERAL_PYB11_TARGET_FLAGS )
endif()

#-------------------------------------------------------------------------------
# BlueOS specific flags
#-------------------------------------------------------------------------------
if (DEFINED ENV{SYS_TYPE})
  if ("$ENV{SYS_TYPE}" STREQUAL "blueos_3_ppc64le_ib_p9")
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
      set(CXX_BLUEOS_FLAGS "-Os")    # Needed to prevent relocation overflow errors during link
      add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:${CXX_BLUEOS_FLAGS}>")
      message("-- Adding ${CXX_BLUEOS_FLAGS} to C++ compile flags")
    endif()
  endif()
endif()
#set(CXX_STRIP_FLAGS "-fdata-sections;-ffunction-sections")
#set(CXX_LINK_STRIP_FLAGS "-Wl,--gc-sections")
#set(CXX_LINK_STRIP_FLAGS "-Wl,-z combreloc")
#add_link_options("${CXX_LINK_STRIP_FLAGS}")

