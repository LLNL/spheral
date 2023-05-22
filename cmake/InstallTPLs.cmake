#-----------------------------------------------------------------------------------
# Define the Third Party Libs to be used here
#-----------------------------------------------------------------------------------

# Do NOT add any TPLs to the clean target
set_directory_properties(PROPERTIES CLEAN_NO_CUSTOM 1)

#-----------------------------------------------------------------------------------
# Submodules
#-----------------------------------------------------------------------------------

# Things all Spheral packages have in their include path
set(SPHERAL_EXTERN_INCLUDES )

# PolyClipper
if (NOT polyclipper_DIR)
  set(polyclipper_DIR "${SPHERAL_ROOT_DIR}/extern/PolyClipper" CACHE PATH "")
endif()
set(polyclipper_INCLUDES "${polyclipper_DIR}/src")

list(APPEND SPHERAL_EXTERN_INCLUDES ${polyclipper_INCLUDES})


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

  list(APPEND SPHERAL_EXTERN_INCLUDES ${PYBIND11_ROOT_DIR}/include)
endif()

#-----------------------------------------------------------------------------------
# Find pre-compiled TPLs
#-----------------------------------------------------------------------------------

# Initialize TPL options
include(${SPHERAL_ROOT_DIR}/cmake/spheral/SpheralHandleTPL.cmake)

# These libs are always needed
Spheral_Handle_TPL(zlib spheral_depends cxx)
Spheral_Handle_TPL(boost spheral_depends cxx)
Spheral_Handle_TPL(eigen spheral_depends cxx)
Spheral_Handle_TPL(qhull spheral_depends cxx)
Spheral_Handle_TPL(silo spheral_depends cxx)

# AXOM PUlls in HDF5 and Conduit for us
#Spheral_Handle_TPL(conduit spheral_depends cxx)
find_package(axom REQUIRED QUIET NO_DEFAULT_PATH PATHS ${axom_DIR}/lib/cmake)
if(axom_FOUND)
  list(APPEND spheral_blt_cxx_depends axom axom::fmt)
  blt_patch_target(NAME axom::fmt TREAT_INCLUDES_AS_SYSTEM On)
  message(STATUS "Found axom: ${axom_DIR} (found version ${axom_VERSION})")
endif()

# Axom imports hdf5 lib but Spheral also requires hdf5_hl
Spheral_Handle_TPL(hdf5 spheral_depends cxx)

# Some libraries are optional
if (ENABLE_ANEOS)
  Spheral_Handle_TPL(aneos spheral_depends cxx)
endif()
if (ENABLE_OPENSUBDIV)
  Spheral_Handle_TPL(opensubdiv spheral_depends cxx)
endif()
if(ENABLE_TIMER)
  Spheral_Handle_TPL(caliper spheral_depends cxx)
endif()

# Only needed when building the python interface of spheral
if(NOT ENABLE_CXXONLY)
  Spheral_Handle_TPL(python spheral_depends cxx)
  #Spheral_Handle_TPL(pyb11generator spheral_depends py)
  #Spheral_Handle_TPL(pybind11 spheral_depends py)
endif()

Spheral_Handle_TPL(polytope spheral_depends cxx)
#Spheral_Handle_TPL(polyclipper spheral_depends cxx)

if (EXISTS ${EXTERNAL_SPHERAL_TPL_CMAKE})
  include(${EXTERNAL_SPHERAL_TPL_CMAKE})
endif()
