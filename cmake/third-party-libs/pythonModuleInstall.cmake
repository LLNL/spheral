set(MODULE_CMAKE_DIR ${PROJECT_SOURCE_DIR}/../cmake/Modules)

macro(Install_PipModules lib_name)
  add_custom_target(
    ${lib_name}
    COMMAND ${PIP_EXE} ${OUT_PROTOCOL_PIP} download --no-binary :all -d ${CACHE_DIR} ${${lib_name}}
    COMMAND ${PIP_EXE} ${OUT_PROTOCOL_PIP} install --upgrade ${${lib_name}} --no-index --find-links ${CACHE_DIR}
    DEPENDS pip-install ${${lib_name}_DEPENDS}
  )
endmacro()

macro(Install_CUSTOM PipModules lib_name)
endmacro()

set(pip-setup-modules
  setuptools
  wheel
)

set(pip-modules
  numpy-stl
  PYB11Generator
  mpi4py
  matplotlib
  decorator
  h5py
  sphinx
  sphinx_rtd_theme
  cython
  sobol
  scipy
  pipreqs
)

set(pip-custom-modules
  gnuplot
)

set(pip-modules_DEPENDS pip-setup-modules)

Install_PipModules(pip-setup-modules)
Install_PipModules(pip-modules)

foreach(lib_name ${pip-custom-modules})
  include(${MODULE_CMAKE_DIR}/${lib_name}Install.cmake)
  list(APPEND spheral_py_depends ${lib_name})
endforeach()

