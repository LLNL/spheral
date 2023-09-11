#-----------------------------------------------------------------------------------
# Define the Third Party Libs to be used here
#-----------------------------------------------------------------------------------

# Do NOT add any TPLs to the clean target
set_directory_properties(PROPERTIES CLEAN_NO_CUSTOM 1)

# Set the location of the <tpl>.cmake files
set(TPL_CMAKE_DIR ${SPHERAL_ROOT_DIR}/cmake/tpl)

#-----------------------------------------------------------------------------------
# Submodules
#-----------------------------------------------------------------------------------

# PolyClipper
if (NOT polyclipper_DIR)
  set(polyclipper_DIR "${SPHERAL_ROOT_DIR}/extern/PolyClipper" CACHE PATH "")
endif()
add_subdirectory(${polyclipper_DIR})
list(APPEND spheral_blt_depends PolyClipperAPI)
install(TARGETS PolyClipperAPI
  EXPORT spheral_cxx-targets
  DESTINATION lib/cmake)
set_target_properties(PolyClipperAPI PROPERTIES EXPORT_NAME spheral::PolyClipperAPI)

if (NOT ENABLE_CXXONLY)
  # Find the appropriate Python
  set(Python3_ROOT_DIR ${python_DIR})
  find_package(Python3 COMPONENTS Interpreter Development)

  # Set the PYB11Generator path
  if (NOT PYB11GENERATOR_ROOT_DIR)
    set(PYB11GENERATOR_ROOT_DIR "${SPHERAL_ROOT_DIR}/extern/PYB11Generator" CACHE PATH "")
  endif()
  include(${PYB11GENERATOR_ROOT_DIR}/cmake/PYB11Generator.cmake)

  # Set the pybind11 path
  if (NOT PYBIND11_ROOT_DIR)
    set(PYBIND11_ROOT_DIR "${PYB11GENERATOR_ROOT_DIR}/extern/pybind11" CACHE PATH "")
  endif()

  list(APPEND SPHERAL_SUBMOD_INCLUDES ${PYBIND11_ROOT_DIR}/include)
endif()

set_property(GLOBAL PROPERTY SPHERAL_SUBMOD_INCLUDES "${SPHERAL_SUBMOD_INCLUDES}")
#-----------------------------------------------------------------------------------
# Find pre-compiled TPLs
#-----------------------------------------------------------------------------------

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
if(NOT ENABLE_CXXONLY)
  list(APPEND SPHERAL_EXTERN_LIBS python)
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
  Spheral_Handle_TPL(${lib} FALSE)
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

if (EXISTS ${EXTERNAL_SPHERAL_TPL_CMAKE})
  include(${EXTERNAL_SPHERAL_TPL_CMAKE})
endif()
