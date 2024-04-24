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
  find_package(Python3 COMPONENTS Interpreter Development)
  set(PYTHON_EXE ${Python3_EXECUTABLE})
  set(SPHERAL_SITE_PACKAGES_PATH "lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages" )
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
  blt_convert_to_system_includes(TARGET PolyClipperAPI)
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
  # BLT Macro for doing this
  blt_convert_to_system_includes(TARGET ${fmt_name})
endif()
# This is a hack to handle transitive issues that come
# from using object libraries with newer version of axom
foreach(_comp ${AXOM_COMPONENTS_ENABLED})
  list(APPEND SPHERAL_BLT_DEPENDS axom::${_comp})
  get_target_property(axom_deps axom::${_comp} INTERFACE_LINK_LIBRARIES)
  list(APPEND SPHERAL_BLT_DEPENDS ${axom_deps})
endforeach()

# TPLs that must be imported
list(APPEND SPHERAL_EXTERN_LIBS boost eigen qhull silo hdf5 polytope)

blt_list_append( TO SPHERAL_EXTERN_LIBS ELEMENTS aneos IF ENABLE_ANEOS)
blt_list_append( TO SPHERAL_EXTERN_LIBS ELEMENTS opensubdiv IF ENABLE_OPENSUBDIV)
blt_list_append( TO SPHERAL_EXTERN_LIBS ELEMENTS caliper IF ENABLE_TIMER)

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

# Copied from serac, needed to bypass generator expression issue during export
set(_props)
if( ${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.13.0" )
  list(APPEND _props INTERFACE_LINK_OPTIONS)
endif()
list(APPEND _props INTERFACE_COMPILE_OPTIONS)
foreach(_target axom axom::openmp)
  if(TARGET ${_target})
    message(STATUS "Removing OpenMP Flags from target[${_target}]")
    foreach(_prop ${_props})
      get_target_property(_flags ${_target} ${_prop})
      if ( _flags )
	string( REPLACE "${OpenMP_CXX_FLAGS}" ""
	  correct_flags "${_flags}" )
	string( REPLACE "${OpenMP_Fortran_FLAGS}" ""
	  correct_flags "${correct_flags}" )
	set_target_properties( ${_target} PROPERTIES ${_prop} "${correct_flags}" )
      endif()
    endforeach()
  endif()
endforeach()
