set_directory_properties(PROPERTIES CLEAN_NO_CUSTOM 1)
set(CACHE_DIR ${PROJECT_SOURCE_DIR}/src/tpl/cache)
set(PATCH_DIR ${PROJECT_SOURCE_DIR}/src/tpl/patch)
set(TPL_CMAKE_DIR ${CMAKE_MODULE_PATH}/tpl)

include(DemoCMake)

set(BUILD_TPL ON CACHE BOOL "")
set(ENABLE_CXXONLY OFF CACHE BOOL "")

set(boost_BUILD ON CACHE BOOL "")
set(eigen_BUILD ON CACHE BOOL "")
set(qhull_BUILD ON CACHE BOOL "")
set(polytope_BUILD ON CACHE BOOL "")
set(hdf5_BUILD ON CACHE BOOL "")
set(silo_BUILD ON CACHE BOOL "")

set(pybind11_BUILD ON CACHE BOOL "")
set(python_BUILD ON CACHE BOOL "")
set(pip_BUILD ON CACHE BOOL "")

Spheral_Handle_TPL(boost spheral_depends)
Spheral_Handle_TPL(eigen spheral_depends)
Spheral_Handle_TPL(qhull spheral_depends)
Spheral_Handle_TPL(hdf5 spheral_depends)
Spheral_Handle_TPL(silo spheral_depends)

if(NOT ENABLE_CXXONLY)
  Spheral_Handle_TPL(python spheral_depends)
  Spheral_Handle_TPL(pip spheral_py_depends)
  include(${TPL_CMAKE_DIR}/pythonModuleInstall.cmake)

  Spheral_Handle_TPL(polytope spheral_depends)
  Spheral_Handle_TPL(pybind11 spheral_depends)
endif()
