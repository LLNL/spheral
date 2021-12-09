include(${SPHERAL_ROOT_DIR}/cmake/tpl/util/Install_Python_distutils_library.cmake)
set(virtualenv_IMPORT virtualenv PARENT_SCOPE)
set(virtualenv_DEPENDS setuptools wheel PARENT_SCOPE)
