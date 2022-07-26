set(${lib_name}_NO_INCLUDES On)
include(${SPHERAL_ROOT_DIR}/cmake/tpl/util/Install_Python_distutils_library.cmake)
set(setuptools_IMPORT setuptools PARENT_SCOPE)
set(setuptools_NO_INCLUDES On)
