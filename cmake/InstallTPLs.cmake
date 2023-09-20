#-----------------------------------------------------------------------------------
# Define the Third Party Libs to be used here
#-----------------------------------------------------------------------------------

# Do NOT add any TPLs to the clean target
set_directory_properties(PROPERTIES CLEAN_NO_CUSTOM 1)

# Set the location of the <tpl>.cmake files
set(TPL_SPHERAL_CMAKE_DIR ${SPHERAL_ROOT_DIR}/cmake/tpl)

#-----------------------------------------------------------------------------------
# Submodules
#-----------------------------------------------------------------------------------

if (NOT ENABLE_CXXONLY)
  # Find the appropriate Python
  set(Python3_ROOT_DIR ${python_DIR})
  find_package(Python3 COMPONENTS Interpreter Development)
  set(PYTHON_EXE ${Python3_EXECUTABLE})
  list(APPEND spheral_blt_depends Python3::Python)

  # Set the PYB11Generator path
  if (NOT PYB11GENERATOR_ROOT_DIR)
    set(PYB11GENERATOR_ROOT_DIR "${SPHERAL_ROOT_DIR}/extern/PYB11Generator" CACHE PATH "")
  endif()
  # Set the pybind11 path
  if (NOT PYBIND11_ROOT_DIR)
    set(PYBIND11_ROOT_DIR "${PYB11GENERATOR_ROOT_DIR}/extern/pybind11" CACHE PATH "")
  endif()
  include(${PYB11GENERATOR_ROOT_DIR}/cmake/PYB11Generator.cmake)
  list(APPEND spheral_blt_depends pybind11_headers)
  install(TARGETS pybind11_headers
    EXPORT spheral_cxx-targets
    DESTINATION lib/cmake)
  set_target_properties(pybind11_headers PROPERTIES EXPORT_NAME spheral::pybind11_headers)
endif()

# This is currently unfilled in spheral
set_property(GLOBAL PROPERTY SPHERAL_SUBMOD_INCLUDES "${SPHERAL_SUBMOD_INCLUDES}")
#-----------------------------------------------------------------------------------
# Find pre-compiled TPLs
#-----------------------------------------------------------------------------------

# PolyClipper
if (NOT polyclipper_DIR)
  set(polyclipper_DIR "${SPHERAL_ROOT_DIR}/extern/PolyClipper" CACHE PATH "")
endif()
# Must set this so PolyClipper doesn't include unnecessary python scripts
set(POLYCLIPPER_MODULE_GEN OFF)
set(POLYCLIPPER_ENABLE_DOCS OFF)
set(POLYCLIPPER_INSTALL_DIR "PolyClipper/include")
add_subdirectory(${polyclipper_DIR})
# Treat includes as system to prevent warnings
blt_patch_target(NAME PolyClipperAPI TREAT_INCLUDES_AS_SYSTEM ON)
list(APPEND spheral_blt_depends PolyClipperAPI)
install(TARGETS PolyClipperAPI
  EXPORT spheral_cxx-targets
  DESTINATION lib/cmake)
set_target_properties(PolyClipperAPI PROPERTIES EXPORT_NAME spheral::PolyClipperAPI)

# TPLs that can use find_package
list(APPEND SPHERAL_EXTERN_PACKAGES axom)

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

# Initialize TPL options
include(${SPHERAL_ROOT_DIR}/cmake/spheral/SpheralHandleTPL.cmake)

foreach(lib ${SPHERAL_EXTERN_PACKAGES})
  # Not sure why but _lib becomes an empty string after find_package so
  # _lib must be added to list before find package
  list(APPEND spheral_blt_depends ${lib})
  find_package(${lib} REQUIRED QUIET NO_DEFAULT_PATH PATHS ${${lib}_DIR}/lib/cmake)
endforeach()
# Add fmt library to external library list
list(APPEND spheral_blt_depends fmt)
blt_patch_target(NAME fmt TREAT_INCLUDES_AS_SYSTEM ON)

# Create target library for each external library
foreach(lib ${SPHERAL_EXTERN_LIBS})
  Spheral_Handle_TPL(${lib} ${TPL_SPHERAL_CMAKE_DIR})
  list(APPEND spheral_blt_depends ${lib})
endforeach()

# Install each TPL target library
foreach(lib ${SPHERAL_EXTERN_LIBS})
  get_target_property(_is_imported ${lib} IMPORTED)
  if(NOT ${_is_imported})
    install(TARGETS ${lib}
      EXPORT spheral_cxx-targets
      DESTINATION lib/cmake)
    set_target_properties(${lib} PROPERTIES EXPORT_NAME spheral::${lib})
  endif()
endforeach()
# Note: spheral_blt_depends is made global after this in SetupSpheral.cmake

# This calls LLNLSpheralInstallTPLs.cmake
if (EXISTS ${EXTERNAL_SPHERAL_TPL_CMAKE})
  include(${EXTERNAL_SPHERAL_TPL_CMAKE})
endif()
