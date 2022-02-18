include(${SPHERAL_ROOT_DIR}/cmake/tpl/util/Install_Python_distutils_library.cmake)

#
# Buildtime Python Module Dependencies
#
set(PYTHONENV ${pip_DIR})

set(pip_build_modules
    setuptools
    decorator==4.4.2
    virtualenv==20.2.2
    pyb11generator==1.0.12)

if (ENABLE_DOCS)
  list(APPEND pip_build_modules
    sphinx==1.8.5
    sphinx_rtd_theme==0.5.0)
endif()

foreach(lib ${pip_build_modules})
  string(REGEX REPLACE "[\=]+[^ ]*" "" lib_name_str ${lib})

  # Initialize our libr}ies with the TPL system. This builds up our PYTHONENV variable...
  Spheral_Handle_TPL(${lib_name_str} spheral_depends)
  list(APPEND PYTHONENV ${${lib_name_str}_DIR})
endforeach()

# Create the final build time python environment string.
string(REPLACE ";" ":" BUILDTIME_PYTHONENV_STR "${PYTHONENV}")


# Run the actual installer for the modules.
foreach(lib ${pip_build_modules})
  Install_Pip_Module(${lib})
endforeach()


#
# Runtime Python Module Dependencies
#
set(pip-runtime-modules
    numpy==1.16.6
    numpy-stl==2.11.2
    matplotlib==2.2.5
    decorator==4.4.2
    h5py==2.10.0
    docutils==0.17.1
    twine==1.15.0
    cython==0.29.21
    sobol==0.9
    scipy==1.2.3
    pipreqs==0.4.10
    importlib_metadata==2.1.1
    )

# Only needed when we're allowing MPI parallelism
if (ENABLE_MPI)
  list(APPEND pip-runtime-modules mpi4py==3.0.3)
endif()

# Generate our requirements.txt file for runtime python dependencies
string(REPLACE ";" "\n" pip_rutime_modules_str "${pip-runtime-modules}")
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/tpl/util/requirements.in"
  "${CMAKE_BINARY_DIR}/scripts/requirements.txt"
  )


#
# Custom Runtime Python Module Dependencies
#
set(pip-custom-modules
    gnuplot
    ats)

foreach(lib_name ${pip-custom-modules})
  include(${TPL_CMAKE_DIR}/${lib_name}.cmake)
endforeach()
