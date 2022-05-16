set(PYBIND11_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${lib_name})
set(PYBIND11_CACHE "${CACHE_DIR}/v2.4.3.tar.gz")
set(PYBIND11_URL "https://github.com/pybind/pybind11/archive/v2.4.3.tar.gz")

set(${lib_name}_libs )
set(PYBIND11_INSTALL_DIR ${${lib_name}_DIR} PARENT_SCOPE)
