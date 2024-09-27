#-----------------------------------------------------------------------------------
# Define the Third Party Libs to be used here
#-----------------------------------------------------------------------------------

# Do NOT add any TPLs to the clean target
set_directory_properties(PROPERTIES CLEAN_NO_CUSTOM 1)

# Set the location of the <tpl>.cmake files
set(TPL_SPHERAL_CMAKE_DIR ${SPHERAL_ROOT_DIR}/cmake/tpl)

# Initialize TPL options
include(${SPHERAL_ROOT_DIR}/cmake/spheral/SpheralHandleTPL.cmake)
include(${SPHERAL_ROOT_DIR}/cmake/spheral/SpheralHandleExt.cmake)

#-----------------------------------------------------------------------------------
# Submodules
#-----------------------------------------------------------------------------------

if (NOT ENABLE_CXXONLY)
  # Find the appropriate Python
  find_package(Python3 COMPONENTS Interpreter Development)
  set(PYTHON_EXE ${Python3_EXECUTABLE})
  set(SPHERAL_SITE_PACKAGES_PATH "lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages" )
  list(APPEND SPHERAL_CXX_DEPENDS Python3::Python)

  # Set the PYB11Generator path
  if (NOT PYB11GENERATOR_ROOT_DIR)
    set(PYB11GENERATOR_ROOT_DIR "${SPHERAL_ROOT_DIR}/extern/PYB11Generator" CACHE PATH "")
  endif()
  # Set the pybind11 path
  if (NOT PYBIND11_ROOT_DIR)
    set(PYBIND11_ROOT_DIR "${PYB11GENERATOR_ROOT_DIR}/extern/pybind11" CACHE PATH "")
  endif()
  include(${PYB11GENERATOR_ROOT_DIR}/cmake/PYB11Generator.cmake)
  list(APPEND SPHERAL_CXX_DEPENDS pybind11_headers)
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

# Any targets that used find package must be added to these lists
set(SPHERAL_FP_TPLS )
set(SPHERAL_FP_DIRS )

# Use find_package to get axom (which brings in fmt) and patch fmt
find_package(axom REQUIRED NO_DEFAULT_PATH PATHS ${axom_DIR}/lib/cmake)
list(APPEND SPHERAL_BLT_DEPENDS axom )
list(APPEND SPHERAL_FP_TPLS axom)
list(APPEND SPHERAL_FP_DIRS ${axom_DIR}/lib/cmake)

# This is a hack to handle transitive issues that come
# from using object libraries with newer version of axom
foreach(_comp ${AXOM_COMPONENTS_ENABLED})
  get_target_property(axom_deps axom::${_comp} INTERFACE_LINK_LIBRARIES)
  # strip cuda out so we have control over when cuda is enabled
  list(REMOVE_DUPLICATES axom_deps)
  list(REMOVE_ITEM axom_deps cuda)
  blt_convert_to_system_includes(TARGET ${axom_deps})
  list(APPEND SPHERAL_BLT_DEPENDS ${axom_deps})
endforeach()

message("-----------------------------------------------------------------------------")
# Use find_package to get adiak
find_package(adiak REQUIRED NO_DEFAULT_PATH PATHS ${adiak_DIR}/lib/cmake/adiak)
if(adiak_FOUND)
  list(APPEND SPHERAL_BLT_DEPENDS adiak::adiak)
  list(APPEND SPHERAL_FP_TPLS adiak)
  list(APPEND SPHERAL_FP_DIRS ${adiak_DIR})
  message("Found Adiak External Package")
endif()
message("-----------------------------------------------------------------------------")
# Use find_package to get caliper
if (ENABLE_TIMER)
  # Save caliper_DIR because it gets overwritten by find_package
  if(NOT CONFIG_CALIPER_DIR)
    # Only save if it does not exists already
    set(CONFIG_CALIPER_DIR "${caliper_DIR}" CACHE PATH "Configuration Caliper directory")
  endif()
  find_package(caliper REQUIRED NO_DEFAULT_PATH PATHS ${caliper_DIR}/share/cmake/caliper)
  if(caliper_FOUND)
    list(APPEND SPHERAL_BLT_DEPENDS caliper)
    list(APPEND SPHERAL_FP_TPLS caliper)
    list(APPEND SPHERAL_FP_DIRS ${caliper_DIR})
    message("Found Caliper External Package")
  endif()
endif()
message("-----------------------------------------------------------------------------")
find_package(RAJA REQUIRED NO_DEFAULT_PATH PATHS ${raja_DIR})
if (RAJA_FOUND) 
  message("Found RAJA External Package.")
endif()
message("-----------------------------------------------------------------------------")
find_package(umpire REQUIRED NO_DEFAULT_PATH PATHS ${umpire_DIR})
if (umpire_FOUND) 
  message("Found umpire External Package.")
endif()
message("-----------------------------------------------------------------------------")

# Chai
if(chai_DIR AND USE_EXTERNAL_CHAI)
  find_package(chai REQUIRED NO_DEFAULT_PATH PATHS ${chai_DIR})
  if (chai_FOUND) 
    message("Found chai External Package.")
  endif()
  list(APPEND SPHERAL_FP_TPLS chai)
  list(APPEND SPHERAL_FP_DIRS ${chai_DIR})
else()
  message("Using chai Submodule.")
  set(chai_DIR "${SPHERAL_ROOT_DIR}/extern/chai")
  set(CHAI_ENABLE_RAJA_PLUGIN On CACHE BOOL "")
  add_subdirectory(${chai_DIR})
endif()

list(APPEND SPHERAL_BLT_DEPENDS chai camp RAJA umpire)
list(APPEND SPHERAL_FP_TPLS RAJA umpire)
list(APPEND SPHERAL_FP_DIRS ${raja_DIR} ${umpire_DIR})
set_property(GLOBAL PROPERTY SPHERAL_FP_TPLS ${SPHERAL_FP_TPLS})
set_property(GLOBAL PROPERTY SPHERAL_FP_DIRS ${SPHERAL_FP_DIRS})

message("-----------------------------------------------------------------------------")

# TPLs that must be imported
list(APPEND SPHERAL_EXTERN_LIBS boost eigen qhull silo hdf5 polytope)

blt_list_append( TO SPHERAL_EXTERN_LIBS ELEMENTS aneos IF ENABLE_ANEOS)
blt_list_append( TO SPHERAL_EXTERN_LIBS ELEMENTS opensubdiv IF ENABLE_OPENSUBDIV)

# Create and install target library for each external library
foreach(lib ${SPHERAL_EXTERN_LIBS})
  if(NOT TARGET ${lib})
    Spheral_Handle_TPL(${lib} ${TPL_SPHERAL_CMAKE_DIR})
    blt_convert_to_system_includes(TARGET ${lib})
  endif()
  list(APPEND SPHERAL_BLT_DEPENDS ${lib})
endforeach()
# Note: SPHERAL_BLT_DEPENDS is made global after this in SetupSpheral.cmake

# This calls LLNLSpheralInstallTPLs.cmake
if (EXISTS ${EXTERNAL_SPHERAL_TPL_CMAKE})
  include(${EXTERNAL_SPHERAL_TPL_CMAKE})
endif()
