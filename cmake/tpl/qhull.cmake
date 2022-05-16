set(QHULL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${lib_name})
set(QHULL_URL "https://github.com/qhull/qhull/archive/2019.1.tar.gz")
set(QHULL_CACHE "${CACHE_DIR}/2019.1.tar.gz")
set(QHULL_SRC_DIR ${QHULL_PREFIX}/src/qhull/src)

# Setting this to just the release library until we support TPL debug builds on LC
set(${lib_name}_libs libqhull.so)
if(ENABLE_STATIC_TPL)
  set(${lib_name}_libs libqhullstatic.a)
endif()
set(QHULL_BUILD_TYPE Release)
# if (CMAKE_BUILD_TYPE STREQUAL "Debug")
#   set(${lib_name}_libs libqhullstatic_d.a)
# else()
#   set(${lib_name}_libs libqhullstatic.a)
# endif()
# set(QHULL_BUILD_TYPE ${CMAKE_BUILD_TYPE})

