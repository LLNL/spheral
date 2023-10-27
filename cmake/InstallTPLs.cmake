#-----------------------------------------------------------------------------------
# Define the Third Party Libs to be used here
#-----------------------------------------------------------------------------------

# Do NOT add any TPLs to the clean target
set_directory_properties(PROPERTIES CLEAN_NO_CUSTOM 1)

# Set the location of the <tpl>.cmake files
set(TPL_SPHERAL_CMAKE_DIR ${SPHERAL_ROOT_DIR}/cmake/tpl)

# Initialize TPL options
include(${SPHERAL_ROOT_DIR}/cmake/spheral/SpheralHandleTPL.cmake)

#-----------------------------------------------------------------------------------
# Submodules
#-----------------------------------------------------------------------------------

if (NOT ENABLE_CXXONLY)
  # Find the appropriate Python
  set(Python3_ROOT_DIR ${python_DIR})
  find_package(Python3 COMPONENTS Interpreter Development)
  set(PYTHON_EXE ${Python3_EXECUTABLE})
  list(APPEND SPHERAL_BLT_DEPENDS Python3::Python)

  # Set the PYB11Generator path
  if (NOT PYB11GENERATOR_ROOT_DIR)
    set(PYB11GENERATOR_ROOT_DIR "${SPHERAL_ROOT_DIR}/extern/PYB11Generator" CACHE PATH "")
  endif()
  # Set the pybind11 path
  if (NOT PYBIND11_ROOT_DIR)
    set(PYBIND11_ROOT_DIR "${PYB11GENERATOR_ROOT_DIR}/extern/pybind11" CACHE PATH "")
  endif()
  include(${PYB11GENERATOR_ROOT_DIR}/cmake/PYB11Generator.cmake)
  list(APPEND SPHERAL_BLT_DEPENDS pybind11_headers)
  install(TARGETS pybind11_headers
    EXPORT spheral_cxx-targets
    DESTINATION lib/cmake)
  set_target_properties(pybind11_headers PROPERTIES EXPORT_NAME spheral::pybind11_headers)
endif()

# This is currently unfilled in spheral
set_property(GLOBAL PROPERTY SPHERAL_SUBMOD_INCLUDES "${SPHERAL_SUBMOD_INCLUDES}")

# PolyClipper
if (NOT polyclipper_DIR)
  # If no PolyClipper is specified, build it as an internal target
  set(polyclipper_DIR "${SPHERAL_ROOT_DIR}/extern/PolyClipper")
  # Must set this so PolyClipper doesn't include unnecessary python scripts
  set(POLYCLIPPER_MODULE_GEN OFF)
  set(POLYCLIPPER_ENABLE_DOCS OFF)
  set(POLYCLIPPER_INSTALL_DIR "PolyClipper/include")
  add_subdirectory(${polyclipper_DIR} ${CMAKE_CURRENT_BINARY_DIR}/PolyClipper)
  # Treat includes as system to prevent warnings
  blt_patch_target(NAME PolyClipperAPI TREAT_INCLUDES_AS_SYSTEM ON)
  list(APPEND SPHERAL_BLT_DEPENDS PolyClipperAPI)
  install(TARGETS PolyClipperAPI
    EXPORT spheral_cxx-targets
    DESTINATION lib/cmake)
  set_target_properties(PolyClipperAPI PROPERTIES EXPORT_NAME spheral::PolyClipperAPI)
else()
  Spheral_Handle_TPL(polyclipper ${TPL_SPHERAL_CMAKE_DIR})
  list(APPEND SPHERAL_BLT_DEPENDS polyclipper)
endif()

#-----------------------------------------------------------------------------------
# Find pre-compiled TPLs
#-----------------------------------------------------------------------------------

# Use find_package to get axom (which brings in fmt) and patch fmt
find_package(axom REQUIRED QUIET NO_DEFAULT_PATH PATHS ${axom_DIR}/lib/cmake)
if(axom_FOUND)
  list(APPEND SPHERAL_BLT_DEPENDS axom)
  # Add fmt library to external library list
  set(fmt_name fmt)
  # Newer Axom versions call fmt target axom::fmt
  if(NOT TARGET fmt)
    set(fmt_name axom::fmt)
  endif()
  list(APPEND SPHERAL_BLT_DEPENDS ${fmt_name})
  blt_patch_target(NAME ${fmt_name} TREAT_INCLUDES_AS_SYSTEM ON)
endif()
# Potential axom dependencies
list(APPEND AXOM_DEPS umpire RAJA conduit::conduit)
foreach(lib ${AXOM_DEPS})
  if(TARGET ${lib})
    list(APPEND SPHERAL_BLT_DEPENDS ${lib})
  endif()
endforeach()

# TPLs that must be imported
list(APPEND SPHERAL_EXTERN_LIBS zlib boost eigen qhull silo hdf5 polytope)

if(ENABLE_ANEOS)
  list(APPEND SPHERAL_EXTERN_LIBS aneos)
endif()
if(ENABLE_OPENSUBDIV)
  list(APPEND SPHERAL_EXTERN_LIBS opensubdiv)
endif()
if(ENABLE_TIMER)
  list(APPEND SPHERAL_EXTERN_LIBS caliper)
endif()

# Create and install target library for each external library
foreach(lib ${SPHERAL_EXTERN_LIBS})
  if(NOT TARGET ${lib})
    Spheral_Handle_TPL(${lib} ${TPL_SPHERAL_CMAKE_DIR})
  endif()
  list(APPEND SPHERAL_BLT_DEPENDS ${lib})
endforeach()
# Note: SPHERAL_BLT_DEPENDS is made global after this in SetupSpheral.cmake

# This calls LLNLSpheralInstallTPLs.cmake
if (EXISTS ${EXTERNAL_SPHERAL_TPL_CMAKE})
  include(${EXTERNAL_SPHERAL_TPL_CMAKE})
endif()
