set(ZLIB_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${lib_name})
set(ZLIB_DIST "zlib-1.2.11.tar.gz")
set(ZLIB_CACHE "${CACHE_DIR}/${ZLIB_DIST}")
set(ZLIB_URL "http://zlib.net/${ZLIB_DIST}")
set(ZLIB_SRC_DIR ${ZLIB_PREFIX}/src/zlib)

set(${lib_name}_libs libz.a)

set(ZLIB_INSTALL_DIR ${${lib_name}_DIR} PARENT_SCOPE)
