#-----------------------------------------------------------------------------------
# Define the list of Third Party Libs to be installed here
#-----------------------------------------------------------------------------------

# Do NOT add any TPLs to the clean target
set_directory_properties(PROPERTIES CLEAN_NO_CUSTOM 1)

# Initialize TPL options
#include(${SPHERAL_ROOT_DIR}/cmake/spheral/SpheralHandleTPL.cmake)
#set(CMAKE_MODULE_PATH  "${SPHERAL_ROOT_DIR}/cmake/tpl/Modules/")

# These libs are always needed
#Spheral_Handle_TPL(zlib spheral_depends)
set(ZLIB_ROOT ${zlib_DIR})
find_package(ZLIB REQUIRED)
list(APPEND spheral_blt_depends ZLIB::ZLIB)

#Spheral_Handle_TPL(boost spheral_depends)
set(BOOST_ROOT ${boost_DIR})
find_package(Boost REQUIRED COMPONENTS filesystem)
list(APPEND spheral_blt_depends Boost::headers)
list(APPEND spheral_blt_depends Boost::filesystem)

#Spheral_Handle_TPL(eigen spheral_depends)
find_package(Eigen3 REQUIRED PATH ${eigen_DIR}/share/eigen3/cmake)
list(APPEND spheral_blt_depends Eigen3::Eigen)
if(EIGEN3_FOUND)
  message(STATUS "Found Eigen3: ${eigen_DIR} (found version ${EIGEN3_VERSION_STRING})")
endif()

#Spheral_Handle_TPL(qhull spheral_depends)
find_package(Qhull REQUIRED PATHS ${qhull_DIR})
list(APPEND spheral_blt_depends Qhull::libqhull)
if(Qhull_FOUND)
  message(STATUS "Found Qhull: ${qhull_DIR} (found version ${Qhull_VERSION})")
endif()

#Spheral_Handle_TPL(hdf5 spheral_depends)
#set(HDF5_USE_STATIC_LIBRARIES On)
#set(HDF5_ROOT ${hdf5_DIR})
#find_package(HDF5 REQUIRED
#             COMPONENTS C HL
#             )
#             #PATHS ${hdf5_DIR}/share/cmake)
#list(APPEND spheral_blt_depends hdf5::hdf5 hdf5::hdf5_hl)
#message("")

#Spheral_Handle_TPL(conduit spheral_depends)

#Spheral_Handle_TPL(axom spheral_depends)
# AXOM PUlls in HDF5 and Conduit for us
find_package(axom REQUIRED QUIET PATHS ${axom_DIR}/lib/cmake)
list(APPEND spheral_blt_depends axom)
if(axom_FOUND)
  message(STATUS "Found axom: ${axom_DIR} (found version ${axom_VERSION})")
endif()


#Spheral_Handle_TPL(silo spheral_depends)
set(Silo_ROOT ${silo_DIR})
find_package(Silo REQUIRED)

# Some libraries are optional
if (ENABLE_ANEOS)
  Spheral_Handle_TPL(aneos spheral_depends)
endif()
if (ENABLE_OPENSUBDIV)
  Spheral_Handle_TPL(opensubdiv spheral_depends)
endif()

# Only needed when building the python interface of spheral
if(NOT ENABLE_CXXONLY)
  Spheral_Handle_TPL(python spheral_depends)
  Spheral_Handle_TPL(pip spheral_py_depends)
  Spheral_Handle_TPL(pybind11 spheral_depends)

  include(${SPHERAL_ROOT_DIR}/cmake/tpl/pythonModule.cmake)

endif()

Spheral_Handle_TPL(polytope spheral_depends)
Spheral_Handle_TPL(polyclipper spheral_depends)

if (EXISTS ${EXTERNAL_SPHERAL_TPL_CMAKE})
  include(${EXTERNAL_SPHERAL_TPL_CMAKE})
endif()
