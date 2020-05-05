set_directory_properties(PROPERTIES CLEAN_NO_CUSTOM 1)
include(${PROJECT_SOURCE_DIR}/../cmake/DemoCMake.cmake)

set(PYTHON_MODULES
  setuptools
  wheel
  numpy-stl
  PYB11Generator
  mpi4py
  matplotlib
  decorator
  h5py
  sphinx
  sphinx_rtd_theme
  twine
  cython
  sobol
  scipy
  pipreqs
  gnuplot
  )

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

Demo_Handle_TPL(boost spheral_depends)
Demo_Handle_TPL(eigen spheral_depends)
Demo_Handle_TPL(qhull spheral_depends)
Demo_Handle_TPL(hdf5 spheral_depends)
Demo_Handle_TPL(silo spheral_depends)

Demo_Handle_TPL(python spheral_depends)

Demo_Handle_TPL(pip spheral_py_depends)
foreach(module ${PYTHON_MODULES})
  Demo_python_lib(${module} spheral_py_depends)
endforeach()

Demo_Handle_TPL(polytope spheral_depends)

if(NOT ENABLE_CXXONLY)
  Demo_Handle_TPL(pybind11 spheral_depends)
endif()
