set(POLYTOPE_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${lib_name})
set(POLYTOPE_CACHE "${CACHE_DIR}/0.6.2.tar.gz")
set(POLYTOPE_URL "https://github.com/pbtoast/polytope/archive/0.6.2.tar.gz")
set(POLYTOPE_SRC_DIR ${POLYTOPE_PREFIX}/src/polytope)

set(${lib_name}_libs libpolytope.a)

if(NOT ENABLE_CXXONLY)
  list(APPEND POLYTOPE_DEPENDS python-install pyb11generator ${spheral_py_depends})
endif()
