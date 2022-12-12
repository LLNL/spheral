#-----------------------------------------------------------------------------------
# Define the list of Third Party Libs to be installed here
#-----------------------------------------------------------------------------------

# Do NOT add any TPLs to the clean target
set_directory_properties(PROPERTIES CLEAN_NO_CUSTOM 1)

# Initialize TPL options
include(${SPHERAL_ROOT_DIR}/cmake/spheral/SpheralHandleTPL.cmake)

# These libs are always needed
Spheral_Handle_TPL(zlib spheral_depends cxx)
Spheral_Handle_TPL(boost spheral_depends cxx)
Spheral_Handle_TPL(eigen spheral_depends cxx)
Spheral_Handle_TPL(qhull spheral_depends cxx)
Spheral_Handle_TPL(silo spheral_depends cxx)

# AXOM PUlls in HDF5 and Conduit for us
#Spheral_Handle_TPL(hdf5 spheral_depends cxx)
#Spheral_Handle_TPL(conduit spheral_depends cxx)
find_package(axom REQUIRED QUIET NO_DEFAULT_PATH PATHS ${axom_DIR}/lib/cmake)
if(axom_FOUND)
  list(APPEND spheral_blt_cxx_depends axom fmt)
  blt_patch_target(NAME fmt TREAT_INCLUDES_AS_SYSTEM On)
  message(STATUS "Found axom: ${axom_DIR} (found version ${axom_VERSION})")
endif()

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
  Spheral_Handle_TPL(python spheral_depends py)
  Spheral_Handle_TPL(pybind11 spheral_depends py)
  list(APPEND spheral_blt_cxx_depends blt_python)
endif()

Spheral_Handle_TPL(polytope spheral_depends cxx)
Spheral_Handle_TPL(polyclipper spheral_depends cxx)

if (EXISTS ${EXTERNAL_SPHERAL_TPL_CMAKE})
  include(${EXTERNAL_SPHERAL_TPL_CMAKE})
endif()
