# Gneral rule for installing pip modules. Downloads necessary tar/wheel files to the cache dir
# then installs in the future from these files so that offline builds are possible.
macro(Install_PipModules lib_name)
  add_custom_target(
    ${lib_name}
    COMMAND ${PIP_EXE} ${OUT_PROTOCOL_PIP} download --no-binary :all -d ${CACHE_DIR} ${${lib_name}}
    COMMAND ${PIP_EXE} ${OUT_PROTOCOL_PIP} install --upgrade ${${lib_name}} --no-index --find-links ${CACHE_DIR}
    DEPENDS pip-install ${${lib_name}_DEPENDS}
  )
endmacro()

# Pip setup modules these need to be 
# installed before everything else
set(pip-setup-modules
    setuptools
    wheel)

# General pip modules, anything from pipy can 
# be added to this list to install
set(pip-modules
    numpy-stl
    PYB11Generator
    matplotlib
    decorator
    h5py
    sphinx
    sphinx_rtd_theme
    twine
    cython
    sobol
    scipy
    pipreqs)

# Only needed when we're allowing MPI parallelism
if (ENABLE_MPI)
  list(APPEND pip-modules mpi4py)
endif()

# Custom python modules to install that need their
# own install process these will be installed last 
# as their own dedicated targets
set(pip-custom-modules
    gnuplot
    ats)

# pip modules need to depend on the pip-setup-modules 
# being installed before the others
set(pip-modules_DEPENDS pip-setup-modules)

# Run install macro for pipy modules
Install_PipModules(pip-setup-modules)
Install_PipModules(pip-modules)

# Run install for custom modules with thier own .cmake file
foreach(lib_name ${pip-custom-modules})
  include(${TPL_CMAKE_DIR}/${lib_name}.cmake)
endforeach()

