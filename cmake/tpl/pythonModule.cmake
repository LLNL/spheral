# Gneral rule for installing pip modules. Downloads necessary tar/wheel files to the cache dir
# then installs in the future from these files so that offline builds are possible.
macro(Install_PipModules lib_list)
  set(${lib_list}_stamp_file "${CACHE_DIR}/.${lib_list}_pip_download.stamp")
  message("-- Adding pip downloads for ${lib_list} : ${${lib_list}_stamp_file}")
  message("------ ${PIP_EXE} ${OUT_PROTOCOL_PIP} download --no-binary :all -d ${CACHE_DIR} ${${lib_list}}")
  add_custom_target(
    ${lib_list}_stamp_file
    COMMAND echo "-- pip downloading ${lib_list}"
    COMMAND ${PIP_EXE} ${OUT_PROTOCOL_PIP} download --no-binary :all -d ${CACHE_DIR} ${${lib_list}}
    COMMAND touch ${${lib_list}_stamp_file}
    DEPENDS pip-install ${${lib_list}_DEPENDS}
  )
  add_custom_target(
    ${lib_list}
    COMMAND echo "-- pip installing ${lib_list}"
    COMMAND ${PIP_EXE} ${OUT_PROTOCOL_PIP} install --upgrade ${${lib_list}} --no-index --find-links ${CACHE_DIR}
    DEPENDS pip-install ${${lib_list}_DEPENDS} ${lib_list}_stamp_file
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
    numpy-stl==2.11.2
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
set(pip-modules_stamp_file_DEPENDS pip-setup-modules)

# Run install macro for pipy modules
Install_PipModules(pip-setup-modules)
Install_PipModules(pip-modules)

# Run install for custom modules with thier own .cmake file
foreach(lib_name ${pip-custom-modules})
  include(${TPL_CMAKE_DIR}/${lib_name}.cmake)
  list(APPEND spheral_py_depends ${lib_name})
endforeach()

