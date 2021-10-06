include(${SPHERAL_ROOT_DIR}/cmake/tpl/util/MD5.cmake)

if (LC_TPL_DIR)
  set(zlib_BUILD OFF CACHE BOOL "Option to build zlib")
  set(boost_BUILD OFF CACHE BOOL "Option to build boost")
  set(eigen_BUILD OFF CACHE BOOL "Option to build eigen")
  set(qhull_BUILD OFF CACHE BOOL "Option to build qhull")
  set(polytope_BUILD OFF CACHE BOOL "Option to build polytope")
  set(hdf5_BUILD OFF CACHE BOOL "Option to build hdf5")
  set(silo_BUILD OFF CACHE BOOL "Option to build silo")
  set(opensubdiv_BUILD OFF CACHE BOOL "Option to build Opensubdiv")
  set(aneos_BUILD OFF CACHE BOOL "Option to build ANEOS third party lib")
  set(conduit_BUILD OFF CACHE BOOL "Option to build Conduit")
  set(axom_BUILD OFF CACHE BOOL "Option to build Axom")
  set(polyclipper_BUILD OFF CACHE BOOL "Option to build PolyClipper")
  #set(raja_BUILD OFF CACHE BOOL "Option to build Raja")

  set(pybind11_BUILD OFF CACHE BOOL "Option to build pybind11")
  set(python_BUILD OFF CACHE BOOL "Option to build python")
  set(pip_BUILD OFF CACHE BOOL "Option to build pip")

  set(aneos_DIR "${LC_TPL_DIR}/aneos/${ANEOS_MD5}" CACHE STRING "")
  set(axom_DIR "${LC_TPL_DIR}/axom/${AXOM_MD5}" CACHE STRING "")
  set(boost_DIR "${LC_TPL_DIR}/boost/${BOOST_MD5}" CACHE STRING "")
  set(conduit_DIR "${LC_TPL_DIR}/conduit/${CONDUIT_MD5}" CACHE STRING "")
  set(eigen_DIR "${LC_TPL_DIR}/eigen/${EIGEN_MD5}" CACHE STRING "")
  set(hdf5_DIR "${LC_TPL_DIR}/hdf5/${HDF5_MD5}" CACHE STRING "")
  set(opensubdiv_DIR "${LC_TPL_DIR}/opensubdiv/${OPENSUBDIV_MD5}" CACHE STRING "")
  set(polytope_DIR "${LC_TPL_DIR}/polytope/${POLYTOPE_MD5}" CACHE STRING "")
  set(pybind11_DIR "${LC_TPL_DIR}/pybind11/${PYBIND11_MD5}" CACHE STRING "")
  set(python_DIR "${LC_TPL_DIR}/python/${PYTHON_MD5}" CACHE STRING "")
  set(qhull_DIR "${LC_TPL_DIR}/qhull/${QHULL_MD5}" CACHE STRING "")
  set(silo_DIR "${LC_TPL_DIR}/silo/${SILO_MD5}" CACHE STRING "")
  set(zlib_DIR "${LC_TPL_DIR}/zlib/${ZLIB_MD5}" CACHE STRING "")
  set(polyclipper_DIR "${LC_TPL_DIR}/polyclipper/${POLYCLIPPER_MD5}" CACHE STRING "")
  #set(raja_DIR "${LC_TPL_DIR}/raja/${RAJA_MD5}" CACHE STRING "")

  #----------------------------------------------------------------------------
  # We need to build some configurations localy until a better 
  # LC TPL management system is implemented.

  # Build axom locally when MPI disabled...
  if(NOT ENABLE_MPI)
    set(axom_BUILD ON)
    set(axom_DIR "")
  endif()

endif()
