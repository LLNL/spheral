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
  set(polyclipper_INCLUDES "${polyclipper_DIR}/src")
else()
  set(polyclipper_INCLUDES "${polyclipper_DIR}/include")
endif()

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
include(CMakeFindDependencyMacro)
find_dependency(axom REQUIRED NO_DEFAULT_PATH PATHS ${axom_DIR}/lib/cmake)
foreach(_tgt umpire RAJA mfem conduit
             axom::c2c axom::fmt axom::core axom::slic
             axom::sidre axom::mint axom::quest axom)
  if(TARGET ${_tgt})
    list(APPEND spheral_blt_cxx_depends ${_tgt})
    # Do we need the following line?
    # blt_patch_target(NAME ${_tgt} TREAT_INCLUDES_AS_SYSTEM On)
    message(WARNING "${_tgt} is a target.")
  else()
    message(WARNING "${_tgt} **IS NOT** a target.")
  endif()
endforeach()

blt_print_target_properties(TARGET umpire)

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
