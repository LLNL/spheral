include(${SPHERAL_ROOT_DIR}/cmake/tpl/util/Install_Python_distutils_library.cmake)
set(sphinx_IMPORT PYB11Generator PARENT_SCOPE)
set(sphinx_DEPENDS sphinx setuptools PARENT_SCOPE)
